//! Copying and clearing clipboard contents.

use std::collections::hash_map::Entry;
use std::collections::{HashMap, HashSet};
use std::ffi::OsString;
use std::fs::{remove_dir_all, File};
use std::io::{self, Read, Seek, SeekFrom, Write};
use std::path::PathBuf;
use std::pin::Pin;
use std::sync::mpsc::sync_channel;
use std::{iter, thread};

use log::trace;
use rustix::fs::{fcntl_setfl, OFlags};
use wayland_client::globals::GlobalListContents;
use wayland_client::protocol::wl_registry::WlRegistry;
use wayland_client::protocol::wl_seat::WlSeat;
use wayland_client::{
    delegate_dispatch, event_created_child, ConnectError, Dispatch, DispatchError, EventQueue,
};

use crate::common::{self, initialize};
use crate::data_control::{
    self, impl_dispatch_device, impl_dispatch_manager, impl_dispatch_offer, impl_dispatch_source,
};
use crate::seat_data::SeatData;
use crate::utils::is_text;

/// The clipboard to operate on.
#[derive(Copy, Clone, Eq, PartialEq, Debug, Hash, PartialOrd, Ord, Default)]
#[cfg_attr(test, derive(proptest_derive::Arbitrary))]
pub enum ClipboardType {
    /// The regular clipboard.
    #[default]
    Regular,
    /// The "primary" clipboard.
    ///
    /// Working with the "primary" clipboard requires the compositor to support ext-data-control,
    /// or wlr-data-control version 2 or above.
    Primary,
    /// Operate on both clipboards at once.
    ///
    /// Useful for atomically setting both clipboards at once. This option requires the "primary"
    /// clipboard to be supported.
    Both,
}

/// MIME type to offer the copied data under.
#[derive(Clone, Eq, PartialEq, Debug, Hash, PartialOrd, Ord)]
#[cfg_attr(test, derive(proptest_derive::Arbitrary))]
pub enum MimeType {
    /// Detect the MIME type automatically from the data.
    #[cfg_attr(test, proptest(skip))]
    Autodetect,
    /// Offer a number of common plain text MIME types.
    Text,
    /// Offer a specific MIME type.
    Specific(String),
}

/// Source for copying.
#[derive(Clone, Eq, PartialEq, Debug, Hash, PartialOrd, Ord)]
#[cfg_attr(test, derive(proptest_derive::Arbitrary))]
pub enum Source {
    /// Copy contents of the standard input.
    #[cfg_attr(test, proptest(skip))]
    StdIn,
    /// Copy the given bytes.
    Bytes(Box<[u8]>),
}

/// Source for copying, with a MIME type.
///
/// Used for [`copy_multi`].
///
/// [`copy_multi`]: fn.copy_multi.html
#[derive(Clone, Eq, PartialEq, Debug, Hash, PartialOrd, Ord)]
pub struct MimeSource {
    pub source: Source,
    pub mime_type: MimeType,
}

/// Seat to operate on.
#[derive(Clone, Eq, PartialEq, Debug, Hash, PartialOrd, Ord, Default)]
pub enum Seat {
    /// Operate on all existing seats at once.
    #[default]
    All,
    /// Operate on a seat with the given name.
    Specific(String),
}

/// Number of paste requests to serve.
#[derive(Copy, Clone, Eq, PartialEq, Debug, Hash, PartialOrd, Ord, Default)]
pub enum ServeRequests {
    /// Serve requests indefinitely.
    #[default]
    Unlimited,
    /// Serve only the given number of requests.
    Only(usize),
}

/// Options and flags that are used to customize the copying.
#[derive(Clone, Eq, PartialEq, Debug, Default, Hash, PartialOrd, Ord)]
pub struct Options {
    /// The clipboard to work with.
    clipboard: ClipboardType,

    /// The seat to work with.
    seat: Seat,

    /// Trim the trailing newline character before copying.
    ///
    /// This flag is only applied for text MIME types.
    trim_newline: bool,

    /// Do not spawn a separate thread for serving copy requests.
    ///
    /// Setting this flag will result in the call to `copy()` **blocking** until all data sources
    /// it creates are destroyed, e.g. until someone else copies something into the clipboard.
    foreground: bool,

    /// Number of paste requests to serve.
    ///
    /// Limiting the number of paste requests to one effectively clears the clipboard after the
    /// first paste. It can be used when copying e.g. sensitive data, like passwords. Note however
    /// that certain apps may have issues pasting when this option is used, in particular XWayland
    /// clients are known to suffer from this.
    serve_requests: ServeRequests,

    /// Omit additional text mime types which are offered by default if at least one text mime type is provided.
    ///
    /// Omits additionally offered `text/plain;charset=utf-8`, `text/plain`, `STRING`, `UTF8_STRING` and
    /// `TEXT` mime types which are offered by default if at least one text mime type is provided.
    omit_additional_text_mime_types: bool,
}

/// A copy operation ready to start serving requests.
pub struct PreparedCopy {
    queue: EventQueue<State>,
    state: State,
    sources: Vec<data_control::Source>,
}

/// Errors that can occur for copying the source data to a temporary file.
#[derive(thiserror::Error, Debug)]
pub enum SourceCreationError {
    #[error("Couldn't create a temporary directory")]
    TempDirCreate(#[source] io::Error),

    #[error("Couldn't create a temporary file")]
    TempFileCreate(#[source] io::Error),

    #[error("Couldn't copy data to the temporary file")]
    DataCopy(#[source] io::Error),

    #[error("Couldn't write to the temporary file")]
    TempFileWrite(#[source] io::Error),

    #[error("Couldn't open the temporary file for newline trimming")]
    TempFileOpen(#[source] io::Error),

    #[error("Couldn't get the temporary file metadata for newline trimming")]
    TempFileMetadata(#[source] io::Error),

    #[error("Couldn't seek the temporary file for newline trimming")]
    TempFileSeek(#[source] io::Error),

    #[error("Couldn't read the last byte of the temporary file for newline trimming")]
    TempFileRead(#[source] io::Error),

    #[error("Couldn't truncate the temporary file for newline trimming")]
    TempFileTruncate(#[source] io::Error),

    #[cfg(feature = "memfd")]
    #[error("Couldn't create the memfd file")]
    MemFdFileCreate(#[source] memfd::Error),

    #[cfg(feature = "memfd")]
    #[error("Couldn't seal the memfd file")]
    MemFdFileSeal(#[source] memfd::Error),

    #[cfg(feature = "memfd")]
    #[error("File is not compatible with memfd")]
    MemFdFileNotCompatible,
}

/// Errors that can occur for copying and clearing the clipboard.
#[derive(thiserror::Error, Debug)]
pub enum Error {
    #[error("There are no seats")]
    NoSeats,

    #[error("Couldn't open the provided Wayland socket")]
    SocketOpenError(#[source] io::Error),

    #[error("Couldn't connect to the Wayland compositor")]
    WaylandConnection(#[source] ConnectError),

    #[error("Wayland compositor communication error")]
    WaylandCommunication(#[source] DispatchError),

    #[error(
        "A required Wayland protocol ({} version {}) is not supported by the compositor",
        name,
        version
    )]
    MissingProtocol { name: &'static str, version: u32 },

    #[error("The compositor does not support primary selection")]
    PrimarySelectionUnsupported,

    #[error("The requested seat was not found")]
    SeatNotFound,

    #[error("Error copying the source into a temporary file")]
    TempCopy(#[source] SourceCreationError),

    #[error("Couldn't remove the temporary file")]
    TempFileRemove(#[source] io::Error),

    #[error("Couldn't remove the temporary directory")]
    TempDirRemove(#[source] io::Error),

    #[error("Error satisfying a paste request")]
    Paste(#[source] DataSourceError),
}

impl From<common::Error> for Error {
    fn from(x: common::Error) -> Self {
        use common::Error::*;

        match x {
            SocketOpenError(err) => Error::SocketOpenError(err),
            WaylandConnection(err) => Error::WaylandConnection(err),
            WaylandCommunication(err) => Error::WaylandCommunication(err.into()),
            MissingProtocol { name, version } => Error::MissingProtocol { name, version },
        }
    }
}

#[derive(thiserror::Error, Debug)]
pub enum DataSourceError {
    #[error("Couldn't open the data file")]
    FileOpen(#[source] io::Error),

    #[error("Couldn't copy the data to the target file descriptor")]
    Copy(#[source] io::Error),
}

struct State {
    common: common::State,
    got_primary_selection: bool,
    // This bool can be set to true when serving a request: either if an error occurs, or if the
    // number of requests to serve was limited and the last request was served.
    should_quit: bool,
    data_sources: HashMap<String, FileSourceInner>,
    serve_requests: ServeRequests,
    // An error that occurred while serving a request, if any.
    error: Option<DataSourceError>,
}

enum FileSourceInner {
    #[cfg(feature = "memfd")]
    InMemory(File),
    TempFile {
        file: File,
        path: PathBuf,
    },
}

impl FileSourceInner {
    fn duplicate(&self) -> Self {
        match self {
            #[cfg(feature = "memfd")]
            Self::InMemory(file) => {
                Self::InMemory(file.try_clone().expect("unable to clone FD we created"))
            }
            Self::TempFile { file, path } => Self::TempFile {
                file: file.try_clone().expect("unable to clone FD we created"),
                path: path.clone(),
            },
        }
    }
}

struct FileSource {
    mime_type: String,
    source: FileSourceInner,
}

delegate_dispatch!(State: [WlSeat: ()] => common::State);

impl AsMut<common::State> for State {
    fn as_mut(&mut self) -> &mut common::State {
        &mut self.common
    }
}

impl Dispatch<WlRegistry, GlobalListContents> for State {
    fn event(
        _state: &mut Self,
        _proxy: &WlRegistry,
        _event: <WlRegistry as wayland_client::Proxy>::Event,
        _data: &GlobalListContents,
        _conn: &wayland_client::Connection,
        _qhandle: &wayland_client::QueueHandle<Self>,
    ) {
    }
}

impl_dispatch_manager!(State);

impl_dispatch_device!(State, WlSeat, |state: &mut Self, event, seat| {
    match event {
        Event::DataOffer { id } => id.destroy(),
        Event::Finished => {
            state.common.seats.get_mut(seat).unwrap().set_device(None);
        }
        Event::PrimarySelection { .. } => {
            state.got_primary_selection = true;
        }
        _ => (),
    }
});

impl_dispatch_offer!(State);

impl_dispatch_source!(State, |state: &mut Self,
                              source: data_control::Source,
                              event| {
    match event {
        Event::Send { mime_type, fd } => {
            // Check if some other source already handled a paste request and indicated that we should
            // quit.
            if state.should_quit {
                source.destroy();
                return;
            }

            // I'm not sure if it's the compositor's responsibility to check that the mime type is
            // valid. Let's check here just in case.
            let data_source = match state.data_sources.get_mut(&mime_type) {
                Some(source) => source,
                None => {
                    return;
                }
            };

            let copy_result = || {
                let source_file = match data_source {
                    #[cfg(feature = "memfd")]
                    FileSourceInner::InMemory(file) => file,
                    FileSourceInner::TempFile { file, .. } => file,
                };
                // Clear O_NONBLOCK, otherwise io::copy() will stop halfway.
                fcntl_setfl(&fd, OFlags::empty())
                    .map_err(io::Error::from)
                    .map_err(DataSourceError::Copy)?;

                let mut target_file = File::from(fd);
                io::copy(source_file, &mut target_file).map_err(DataSourceError::Copy)
            };

            if let Err(err) = copy_result() {
                state.error = Some(err);
            }

            let done = if let ServeRequests::Only(left) = state.serve_requests {
                let left = left.checked_sub(1).unwrap();
                state.serve_requests = ServeRequests::Only(left);
                left == 0
            } else {
                false
            };

            if done || state.error.is_some() {
                state.should_quit = true;
                source.destroy();
            }
        }
        Event::Cancelled => source.destroy(),
        _ => (),
    }
});

impl Options {
    /// Creates a blank new set of options ready for configuration.
    #[inline]
    pub fn new() -> Self {
        Self::default()
    }

    /// Sets the clipboard to work with.
    #[inline]
    pub fn clipboard(&mut self, clipboard: ClipboardType) -> &mut Self {
        self.clipboard = clipboard;
        self
    }

    /// Sets the seat to use for copying.
    #[inline]
    pub fn seat(&mut self, seat: Seat) -> &mut Self {
        self.seat = seat;
        self
    }

    /// Sets the flag for trimming the trailing newline.
    ///
    /// This flag is only applied for text MIME types.
    #[inline]
    pub fn trim_newline(&mut self, trim_newline: bool) -> &mut Self {
        self.trim_newline = trim_newline;
        self
    }

    /// Sets the flag for not spawning a separate thread for serving copy requests.
    ///
    /// Setting this flag will result in the call to `copy()` **blocking** until all data sources
    /// it creates are destroyed, e.g. until someone else copies something into the clipboard.
    #[inline]
    pub fn foreground(&mut self, foreground: bool) -> &mut Self {
        self.foreground = foreground;
        self
    }

    /// Sets the number of requests to serve.
    ///
    /// Limiting the number of requests to one effectively clears the clipboard after the first
    /// paste. It can be used when copying e.g. sensitive data, like passwords. Note however that
    /// certain apps may have issues pasting when this option is used, in particular XWayland
    /// clients are known to suffer from this.
    #[inline]
    pub fn serve_requests(&mut self, serve_requests: ServeRequests) -> &mut Self {
        self.serve_requests = serve_requests;
        self
    }

    /// Sets the flag for omitting additional text mime types which are offered by default if at least one text mime type is provided.
    ///
    /// Omits additionally offered `text/plain;charset=utf-8`, `text/plain`, `STRING`, `UTF8_STRING` and
    /// `TEXT` mime types which are offered by default if at least one text mime type is provided.
    #[inline]
    pub fn omit_additional_text_mime_types(
        &mut self,
        omit_additional_text_mime_types: bool,
    ) -> &mut Self {
        self.omit_additional_text_mime_types = omit_additional_text_mime_types;
        self
    }

    /// Invokes the copy operation. See `copy()`.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # extern crate wl_clipboard_rs;
    /// # use wl_clipboard_rs::copy::Error;
    /// # fn foo() -> Result<(), Error> {
    /// use wl_clipboard_rs::copy::{MimeType, Options, Source};
    ///
    /// let opts = Options::new();
    /// opts.copy(Source::Bytes([1, 2, 3][..].into()), MimeType::Autodetect)?;
    /// # Ok(())
    /// # }
    /// ```
    #[inline]
    pub fn copy(self, source: Source, mime_type: MimeType) -> Result<(), Error> {
        copy(self, source, mime_type)
    }

    /// Invokes the copy_multi operation. See `copy_multi()`.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # extern crate wl_clipboard_rs;
    /// # use wl_clipboard_rs::copy::Error;
    /// # fn foo() -> Result<(), Error> {
    /// use wl_clipboard_rs::copy::{MimeSource, MimeType, Options, Source};
    ///
    /// let opts = Options::new();
    /// opts.copy_multi(vec![MimeSource { source: Source::Bytes([1, 2, 3][..].into()),
    ///                                   mime_type: MimeType::Autodetect },
    ///                      MimeSource { source: Source::Bytes([7, 8, 9][..].into()),
    ///                                   mime_type: MimeType::Text }])?;
    /// # Ok(())
    /// # }
    /// ```
    #[inline]
    pub fn copy_multi(self, sources: Vec<MimeSource>) -> Result<(), Error> {
        copy_multi(self, sources)
    }

    /// Invokes the prepare_copy operation. See `prepare_copy()`.
    ///
    /// # Panics
    ///
    /// Panics if `foreground` is `false`.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # extern crate wl_clipboard_rs;
    /// # use wl_clipboard_rs::copy::Error;
    /// # fn foo() -> Result<(), Error> {
    /// use wl_clipboard_rs::copy::{MimeSource, MimeType, Options, Source};
    ///
    /// let mut opts = Options::new();
    /// opts.foreground(true);
    /// let prepared_copy = opts.prepare_copy(Source::Bytes([1, 2, 3][..].into()),
    ///                                       MimeType::Autodetect)?;
    /// prepared_copy.serve()?;
    ///
    /// # Ok(())
    /// # }
    /// ```
    #[inline]
    pub fn prepare_copy(self, source: Source, mime_type: MimeType) -> Result<PreparedCopy, Error> {
        prepare_copy(self, source, mime_type)
    }

    /// Invokes the prepare_copy_multi operation. See `prepare_copy_multi()`.
    ///
    /// # Panics
    ///
    /// Panics if `foreground` is `false`.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # extern crate wl_clipboard_rs;
    /// # use wl_clipboard_rs::copy::Error;
    /// # fn foo() -> Result<(), Error> {
    /// use wl_clipboard_rs::copy::{MimeSource, MimeType, Options, Source};
    ///
    /// let mut opts = Options::new();
    /// opts.foreground(true);
    /// let prepared_copy =
    ///     opts.prepare_copy_multi(vec![MimeSource { source: Source::Bytes([1, 2, 3][..].into()),
    ///                                               mime_type: MimeType::Autodetect },
    ///                                  MimeSource { source: Source::Bytes([7, 8, 9][..].into()),
    ///                                               mime_type: MimeType::Text }])?;
    /// prepared_copy.serve()?;
    ///
    /// # Ok(())
    /// # }
    /// ```
    #[inline]
    pub fn prepare_copy_multi(self, sources: Vec<MimeSource>) -> Result<PreparedCopy, Error> {
        prepare_copy_multi(self, sources)
    }
}

impl PreparedCopy {
    /// Starts serving copy requests.
    ///
    /// This function **blocks** until all requests are served or the clipboard is taken over by
    /// some other application.
    pub fn serve(mut self) -> Result<(), Error> {
        // Loop until we're done.
        while !self.state.should_quit {
            self.queue
                .blocking_dispatch(&mut self.state)
                .map_err(Error::WaylandCommunication)?;

            // Check if all sources have been destroyed.
            let all_destroyed = self.sources.iter().all(|x| !x.is_alive());
            if all_destroyed {
                self.state.should_quit = true;
            }
        }

        // Clean up the temp file and directory.
        //
        // We want to try cleaning up all files and folders, so if any errors occur in process,
        // collect them into a vector without interruption, and then return the first one.
        let mut results = Vec::new();
        let mut disk_source_dropped = HashSet::new();
        for data_source in self.state.data_sources.into_values() {
            match data_source {
                // memfd is reference counted, there's no chance of closing something twice.
                #[cfg(feature = "memfd")]
                FileSourceInner::InMemory(file) => {
                    drop(file);
                    results.push(Ok(()));
                }
                FileSourceInner::TempFile { file, mut path } => {
                    drop(file);
                    path.pop();
                    // `data_sources` can contain the same source for many mimes, but we want to free each once.
                    // As every temporary file is under a different `$TMPDIR/` subdirectory, we are
                    // able to dedup based of the subdirectory instead.
                    if disk_source_dropped.contains(&path) {
                        continue;
                    };

                    results.push(remove_dir_all(&path).map_err(Error::TempDirRemove));
                    disk_source_dropped.insert(path);
                }
            }
        }

        // Return the error, if any.
        let result: Result<(), _> = results.into_iter().collect();
        result?;

        // Check if an error occurred during data transfer.
        if let Some(err) = self.state.error.take() {
            return Err(Error::Paste(err));
        }

        Ok(())
    }
}

#[cfg(feature = "memfd")]
fn make_source(
    source: Source,
    mime_type: MimeType,
    trim_newline: bool,
) -> Result<FileSource, SourceCreationError> {
    use std::os::fd::{FromRawFd, IntoRawFd};

    use memfd::{FileSeal, Memfd, MemfdOptions};

    let memfd_file = MemfdOptions::default()
        .allow_sealing(true)
        .close_on_exec(true)
        .create("wl-clipboard")
        .map_err(SourceCreationError::MemFdFileCreate)?;

    let mut memfd_file = memfd_file.into_file();
    let mime_type = get_file_mime_type(&memfd_file, mime_type)?;
    let pinned_file = PinnedFile {
        inner: &mut memfd_file,
        _pinned: std::marker::PhantomPinned,
    };
    let pinned_file = std::pin::pin!(pinned_file);

    write_source_to_file(source, trim_newline, pinned_file, &mime_type)?;

    // SAFETY: Created by `memfd` above and we know nothing in-between
    // swapped out the underlying file descriptor.
    #[allow(unsafe_code)]
    let memfd_file = unsafe { Memfd::from_raw_fd(memfd_file.into_raw_fd()) };

    memfd_file
        .add_seals(&[
            FileSeal::SealWrite,
            FileSeal::SealGrow,
            FileSeal::SealShrink,
            FileSeal::SealSeal,
        ])
        .map_err(SourceCreationError::MemFdFileSeal)?;

    Ok(FileSource {
        mime_type,
        source: FileSourceInner::InMemory(memfd_file.into_file()),
    })
}

#[cfg(not(feature = "memfd"))]
fn make_source(
    source: Source,
    mime_type: MimeType,
    trim_newline: bool,
) -> Result<FileSource, SourceCreationError> {
    let temp_dir = tempfile::tempdir().map_err(SourceCreationError::TempDirCreate)?;
    let mut temp_filename = temp_dir.into_path();
    temp_filename.push("stdin");
    trace!("Temp filename: {}", temp_filename.to_string_lossy());
    // Open the resulting FD as R+W so that both writing the source, and
    // later copying it out to a destination, work as intended.
    let mut temp_file = std::fs::OpenOptions::new()
        .write(true)
        .read(true)
        .create(true)
        .truncate(true)
        .open(&temp_filename)
        .map_err(SourceCreationError::TempFileCreate)?;
    let mime_type = get_file_mime_type(&temp_file, mime_type)?;

    let pinned_file = PinnedFile {
        inner: &mut temp_file,
        _pinned: std::marker::PhantomPinned,
    };
    let pinned_file = std::pin::pin!(pinned_file);

    write_source_to_file(source, trim_newline, pinned_file, &mime_type)?;
    Ok(FileSource {
        mime_type,
        source: FileSourceInner::TempFile {
            file: temp_file,
            path: temp_filename,
        },
    })
}

fn get_file_mime_type(
    source_file: &File,
    mime_type: MimeType,
) -> Result<String, SourceCreationError> {
    let mime_type_str = match mime_type {
        MimeType::Autodetect => match tree_magic_mini::from_file(source_file) {
            Some(magic) => Ok(magic),
            None => Err(SourceCreationError::TempFileOpen(std::io::Error::new(
                std::io::ErrorKind::Other,
                "problem with temp file",
            ))),
        }?
        .to_string(),
        MimeType::Text => "text/plain".to_string(),
        MimeType::Specific(mime_type) => mime_type,
    };
    trace!("Base MIME type: {}", mime_type_str);
    Ok(mime_type_str)
}

// A wrapper which helps force proper `unsafe` annotations
// to ensure FDs don't get swapped out by `write_source_to_file`
// in order to soundly re-consume the FD after.
struct PinnedFile<'a> {
    inner: &'a mut File,
    _pinned: std::marker::PhantomPinned,
}

fn write_source_to_file(
    source: Source,
    trim_newline: bool,
    output_place: Pin<&mut PinnedFile>,
    mime_type: &str,
) -> Result<(), SourceCreationError> {
    // SAFETY: We only write to the file, we never replace it.
    #[allow(unsafe_code)]
    let output_place = unsafe { &mut *output_place.get_unchecked_mut().inner };

    if let Source::Bytes(data) = source {
        output_place
            .write_all(&data)
            .map_err(SourceCreationError::TempFileWrite)?;
    } else {
        // Copy the standard input into the target file.
        io::copy(&mut io::stdin(), output_place).map_err(SourceCreationError::DataCopy)?;
    }

    // Trim the trailing newline if needed.
    if trim_newline && is_text(mime_type) {
        let metadata = output_place
            .metadata()
            .map_err(SourceCreationError::TempFileMetadata)?;
        let length = metadata.len();
        if length > 0 {
            output_place
                .seek(SeekFrom::End(-1))
                .map_err(SourceCreationError::TempFileSeek)?;

            let mut buf = [0];
            output_place
                .read_exact(&mut buf)
                .map_err(SourceCreationError::TempFileRead)?;
            if buf[0] == b'\n' {
                output_place
                    .set_len(length - 1)
                    .map_err(SourceCreationError::TempFileTruncate)?;
            }
        }
    }

    // Restore its offset to the start so that
    // future reads see the previously-written data.
    output_place
        .seek(SeekFrom::Start(0))
        .map(drop)
        .map_err(SourceCreationError::TempFileSeek)
}

fn get_devices(
    primary: bool,
    seat: Seat,
    socket_name: Option<OsString>,
) -> Result<(EventQueue<State>, State, Vec<data_control::Device>), Error> {
    let (mut queue, mut common) = initialize(primary, socket_name)?;

    // Check if there are no seats.
    if common.seats.is_empty() {
        return Err(Error::NoSeats);
    }

    // Go through the seats and get their data devices.
    for (seat, data) in &mut common.seats {
        let device = common
            .clipboard_manager
            .get_data_device(seat, &queue.handle(), seat.clone());
        data.set_device(Some(device));
    }

    let mut state = State {
        common,
        got_primary_selection: false,
        should_quit: false,
        data_sources: HashMap::new(),
        serve_requests: ServeRequests::default(),
        error: None,
    };

    // Retrieve all seat names.
    queue
        .roundtrip(&mut state)
        .map_err(Error::WaylandCommunication)?;

    // Check if the compositor supports primary selection.
    if primary && !state.got_primary_selection {
        return Err(Error::PrimarySelectionUnsupported);
    }

    // Figure out which devices we're interested in.
    let devices = state
        .common
        .seats
        .values()
        .filter_map(|data| {
            let SeatData { name, device, .. } = data;

            let device = device.clone();

            match seat {
                Seat::All => {
                    // If no seat was specified, handle all of them.
                    return device;
                }
                Seat::Specific(ref desired_name) => {
                    if name.as_deref() == Some(desired_name) {
                        return device;
                    }
                }
            }

            None
        })
        .collect::<Vec<_>>();

    // If we didn't find the seat, print an error message and exit.
    //
    // This also triggers when we found the seat but it had no data device; is this what we want?
    if devices.is_empty() {
        return Err(Error::SeatNotFound);
    }

    Ok((queue, state, devices))
}

/// Clears the clipboard for the given seat.
///
/// If `seat` is `None`, clears clipboards of all existing seats.
///
/// # Examples
///
/// ```no_run
/// # extern crate wl_clipboard_rs;
/// # use wl_clipboard_rs::copy::Error;
/// # fn foo() -> Result<(), Error> {
/// use wl_clipboard_rs::{copy::{clear, ClipboardType, Seat}};
///
/// clear(ClipboardType::Regular, Seat::All)?;
/// # Ok(())
/// # }
/// ```
#[inline]
pub fn clear(clipboard: ClipboardType, seat: Seat) -> Result<(), Error> {
    clear_internal(clipboard, seat, None)
}

pub(crate) fn clear_internal(
    clipboard: ClipboardType,
    seat: Seat,
    socket_name: Option<OsString>,
) -> Result<(), Error> {
    let primary = clipboard != ClipboardType::Regular;
    let (mut queue, mut state, devices) = get_devices(primary, seat, socket_name)?;

    for device in devices {
        if clipboard == ClipboardType::Primary || clipboard == ClipboardType::Both {
            device.set_primary_selection(None);
        }
        if clipboard == ClipboardType::Regular || clipboard == ClipboardType::Both {
            device.set_selection(None);
        }
    }

    // We're clearing the clipboard so just do one roundtrip and quit.
    queue
        .roundtrip(&mut state)
        .map_err(Error::WaylandCommunication)?;

    Ok(())
}

/// Prepares a data copy to the clipboard.
///
/// The data is copied from `source` and offered in the `mime_type` MIME type. See `Options` for
/// customizing the behavior of this operation.
///
/// This function can be used instead of `copy()` when it's desirable to separately prepare the
/// copy operation, handle any errors that this may produce, and then start the serving loop,
/// potentially past a fork (which is how `wl-copy` uses it). It is meant to be used in the
/// foreground mode and does not spawn any threads.
///
/// # Panics
///
/// Panics if `foreground` is `false`.
///
/// # Examples
///
/// ```no_run
/// # extern crate wl_clipboard_rs;
/// # use wl_clipboard_rs::copy::Error;
/// # fn foo() -> Result<(), Error> {
/// use wl_clipboard_rs::copy::{MimeSource, MimeType, Options, Source};
///
/// let mut opts = Options::new();
/// opts.foreground(true);
/// let prepared_copy = opts.prepare_copy(Source::Bytes([1, 2, 3][..].into()),
///                                       MimeType::Autodetect)?;
/// prepared_copy.serve()?;
///
/// # Ok(())
/// # }
/// ```
#[inline]
pub fn prepare_copy(
    options: Options,
    source: Source,
    mime_type: MimeType,
) -> Result<PreparedCopy, Error> {
    assert!(options.foreground);

    let sources = vec![MimeSource { source, mime_type }];

    prepare_copy_internal(options, sources, None)
}

/// Prepares a data copy to the clipboard, offering multiple data sources.
///
/// The data from each source in `sources` is copied and offered in the corresponding MIME type.
/// See `Options` for customizing the behavior of this operation.
///
/// If multiple sources specify the same MIME type, the first one is offered. If one of the MIME
/// types is text, all automatically added plain text offers will fall back to the first source
/// with a text MIME type.
///
/// This function can be used instead of `copy()` when it's desirable to separately prepare the
/// copy operation, handle any errors that this may produce, and then start the serving loop,
/// potentially past a fork (which is how `wl-copy` uses it). It is meant to be used in the
/// foreground mode and does not spawn any threads.
///
/// # Panics
///
/// Panics if `foreground` is `false`.
///
/// # Examples
///
/// ```no_run
/// # extern crate wl_clipboard_rs;
/// # use wl_clipboard_rs::copy::Error;
/// # fn foo() -> Result<(), Error> {
/// use wl_clipboard_rs::copy::{MimeSource, MimeType, Options, Source};
///
/// let mut opts = Options::new();
/// opts.foreground(true);
/// let prepared_copy =
///     opts.prepare_copy_multi(vec![MimeSource { source: Source::Bytes([1, 2, 3][..].into()),
///                                               mime_type: MimeType::Autodetect },
///                                  MimeSource { source: Source::Bytes([7, 8, 9][..].into()),
///                                               mime_type: MimeType::Text }])?;
/// prepared_copy.serve()?;
///
/// # Ok(())
/// # }
/// ```
#[inline]
pub fn prepare_copy_multi(
    options: Options,
    sources: Vec<MimeSource>,
) -> Result<PreparedCopy, Error> {
    assert!(options.foreground);

    prepare_copy_internal(options, sources, None)
}

fn prepare_copy_internal(
    options: Options,
    sources: Vec<MimeSource>,
    socket_name: Option<OsString>,
) -> Result<PreparedCopy, Error> {
    let Options {
        clipboard,
        seat,
        trim_newline,
        serve_requests,
        ..
    } = options;

    let primary = clipboard != ClipboardType::Regular;
    let (queue, mut state, devices) = get_devices(primary, seat, socket_name)?;

    state.serve_requests = serve_requests;

    // Collect the source data to copy.
    state.data_sources = {
        let mut data_sources = HashMap::new();
        let mut text_data = None;

        for MimeSource { source, mime_type } in sources.into_iter() {
            let FileSource { mime_type, source } =
                make_source(source, mime_type, trim_newline).map_err(Error::TempCopy)?;

            match data_sources.entry(mime_type) {
                Entry::Occupied(_) => {
                    // This MIME type has already been specified, so ignore it.
                    match source {
                        #[cfg(feature = "memfd")]
                        FileSourceInner::InMemory(file) => drop(file),
                        FileSourceInner::TempFile { file, mut path, .. } => {
                            // Be explicit during cleanup.
                            drop(file);
                            // Get the tempdir parent created prior.
                            path.pop();
                            remove_dir_all(path).map_err(Error::TempFileRemove)?;
                        }
                    }
                }
                Entry::Vacant(entry) => {
                    if !options.omit_additional_text_mime_types
                        && text_data.is_none()
                        && is_text(entry.key())
                    {
                        text_data = Some((entry.into_key(), source));
                    } else {
                        entry.insert(source);
                    }
                }
            }
        }

        // If the MIME type is text, offer it in some other common formats.
        if let Some((original_text_mime, text_data)) = text_data {
            let text_mimes = [
                "text/plain;charset=utf-8",
                "text/plain",
                "STRING",
                "UTF8_STRING",
                "TEXT",
            ];

            // XXX: The code below optimizes for minimal fd duplication/cloning since
            // they are more likely to become exhausted or cause issues then small memory allocations.

            // We don't want to overwrite an explicit mime type, because it might be bound to a
            // different `data_source`.
            let sources_needed = text_mimes
                .into_iter()
                .chain(Some(original_text_mime.as_str()))
                .filter(|mime| !data_sources.contains_key(*mime))
                .collect::<HashSet<_>>();

            let sources_count = sources_needed.len();
            let mut sources = Vec::with_capacity(sources_count);
            sources
                .extend(std::iter::repeat_with(|| text_data.duplicate()).take(sources_count - 1));
            sources.push(text_data);

            for (mime_type, text_data) in sources_needed.into_iter().zip(sources) {
                data_sources.insert(mime_type.to_string(), text_data);
            }
        }
        data_sources
    };

    // Create an iterator over (device, primary) for source creation later.
    //
    // This is needed because for ClipboardType::Both each device needs to appear twice because
    // separate data sources need to be made for the regular and the primary clipboards (data
    // sources cannot be reused).
    let devices_iter = devices.iter().flat_map(|device| {
        let first = match clipboard {
            ClipboardType::Regular => iter::once((device, false)),
            ClipboardType::Primary => iter::once((device, true)),
            ClipboardType::Both => iter::once((device, false)),
        };

        let second = if clipboard == ClipboardType::Both {
            iter::once(Some((device, true)))
        } else {
            iter::once(None)
        };

        first.chain(second.flatten())
    });

    // Create the data sources and set them as selections.
    let sources = devices_iter
        .map(|(device, primary)| {
            let data_source = state
                .common
                .clipboard_manager
                .create_data_source(&queue.handle());

            for mime_type in state.data_sources.keys() {
                data_source.offer(mime_type.clone());
            }

            if primary {
                device.set_primary_selection(Some(&data_source));
            } else {
                device.set_selection(Some(&data_source));
            }

            // If we need to serve 0 requests, kill the data source right away.
            if let ServeRequests::Only(0) = state.serve_requests {
                data_source.destroy();
            }
            data_source
        })
        .collect::<Vec<_>>();

    Ok(PreparedCopy {
        queue,
        state,
        sources,
    })
}

/// Copies data to the clipboard.
///
/// The data is copied from `source` and offered in the `mime_type` MIME type. See `Options` for
/// customizing the behavior of this operation.
///
/// # Examples
///
/// ```no_run
/// # extern crate wl_clipboard_rs;
/// # use wl_clipboard_rs::copy::Error;
/// # fn foo() -> Result<(), Error> {
/// use wl_clipboard_rs::copy::{copy, MimeType, Options, Source};
///
/// let opts = Options::new();
/// copy(opts, Source::Bytes([1, 2, 3][..].into()), MimeType::Autodetect)?;
/// # Ok(())
/// # }
/// ```
#[inline]
pub fn copy(options: Options, source: Source, mime_type: MimeType) -> Result<(), Error> {
    let sources = vec![MimeSource { source, mime_type }];
    copy_internal(options, sources, None)
}

/// Copies data to the clipboard, offering multiple data sources.
///
/// The data from each source in `sources` is copied and offered in the corresponding MIME type.
/// See `Options` for customizing the behavior of this operation.
///
/// If multiple sources specify the same MIME type, the first one is offered. If one of the MIME
/// types is text, all automatically added plain text offers will fall back to the first source
/// with a text MIME type.
///
/// # Examples
///
/// ```no_run
/// # extern crate wl_clipboard_rs;
/// # use wl_clipboard_rs::copy::Error;
/// # fn foo() -> Result<(), Error> {
/// use wl_clipboard_rs::copy::{MimeSource, MimeType, Options, Source};
///
/// let opts = Options::new();
/// opts.copy_multi(vec![MimeSource { source: Source::Bytes([1, 2, 3][..].into()),
///                                   mime_type: MimeType::Autodetect },
///                      MimeSource { source: Source::Bytes([7, 8, 9][..].into()),
///                                   mime_type: MimeType::Text }])?;
/// # Ok(())
/// # }
/// ```
#[inline]
pub fn copy_multi(options: Options, sources: Vec<MimeSource>) -> Result<(), Error> {
    copy_internal(options, sources, None)
}

pub(crate) fn copy_internal(
    options: Options,
    sources: Vec<MimeSource>,
    socket_name: Option<OsString>,
) -> Result<(), Error> {
    if options.foreground {
        prepare_copy_internal(options, sources, socket_name)?.serve()
    } else {
        // The copy must be prepared on the thread because PreparedCopy isn't Send.
        // To receive errors from prepare_copy, use a channel.
        let (tx, rx) = sync_channel(1);

        thread::spawn(
            move || match prepare_copy_internal(options, sources, socket_name) {
                Ok(prepared_copy) => {
                    // prepare_copy completed successfully, report that.
                    drop(tx.send(None));

                    // There's nobody listening for errors at this point, just drop it.
                    drop(prepared_copy.serve());
                }
                Err(err) => drop(tx.send(Some(err))),
            },
        );

        if let Some(err) = rx.recv().unwrap() {
            return Err(err);
        }

        Ok(())
    }
}
