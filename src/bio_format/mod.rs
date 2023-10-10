pub use nu_protocol::{Span, Value};
/// SAM + BAM parsing facility.
pub mod bam;
/// BCF + VCF parsing facility.
pub mod bcf;
/// BED parsing facility
pub mod bed;
/// CRAM parsing facility.
pub mod cram;
/// Fasta parsing facility.
pub mod fasta;
/// GFA parsing utility
pub mod gfa;
/// GFF(3) parsing facility
pub mod gff;

/// Compression enum
pub enum Compression {
    Uncompressed,
    Gzipped,
}

pub trait SpanExt {
    fn with_string<S: ToString>(&self, s: S) -> Value;
    fn with_string_or<S: ToString>(&self, s: Option<S>, default: &str) -> Value;
}

impl SpanExt for Span {
    fn with_string<S: ToString>(&self, s: S) -> Value {
        Value::String {
            val: s.to_string(),
            span: *self,
        }
    }

    fn with_string_or<S: ToString>(&self, s: Option<S>, default: &str) -> Value {
        Value::String {
            val: s.map(|s| s.to_string()).unwrap_or(default.into()),
            span: *self,
        }
    }
}
