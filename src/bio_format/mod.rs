/// SAM + BAM parsing facility.
pub mod bam;
/// BCF + VCF parsing facility.
pub mod bcf;
/// CRAM parsing facility.
pub mod cram;
/// Fasta parsing facility.
pub mod fasta;
/// GFF(3) parsing facility
pub mod gff;

/// Compression enum
pub enum Compression {
    Uncompressed,
    Gzipped,
}
