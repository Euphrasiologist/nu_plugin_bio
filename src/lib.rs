//! A bioinformatics parsing library for nushell.

/// Where the core parsers live.
mod bio;
/// Handle all the file types.
mod bio_format;
/// Nushell logic handling.
mod nu;

/// Expose the [`Bio`] struct here.
pub use bio::Bio;
