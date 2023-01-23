/// SAM + BAM parsing facility.
pub mod bam;
/// BCF + VCF parsing facility.
pub mod bcf;
/// CRAM parsing facility.
pub mod cram;
/// Fasta parsing facility.
pub mod fasta;
/// GFA parsing utility
pub mod gfa;
/// GFF(3) parsing facility
pub mod gff;
/// BED parsing facility
pub mod bed;

/// Compression enum
pub enum Compression {
    Uncompressed,
    Gzipped,
}
