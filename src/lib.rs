//! A bioinformatics parsing library for nushell.

/// Where the core parsers live.
mod bio;
/// Nushell logic handling.
mod nu;

/// Expose the [`Bio`] struct here.
pub use bio::Bio;
