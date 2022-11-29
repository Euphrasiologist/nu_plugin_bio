use crate::bio_format::Compression;
use crate::Bio;
use nu_plugin::{EvaluatedCall, LabeledError, Plugin};
use nu_protocol::{Category, Signature, Value};

impl Plugin for Bio {
    fn signature(&self) -> Vec<Signature> {
        vec![
            Signature::build("from fasta")
                .usage("Parse a fasta file.")
                .switch(
                    "description",
                    "parse the fasta header description",
                    Some('d'),
                )
                .category(Category::Experimental),
            Signature::build("from fasta.gz")
                .usage("Parse a gzipped fasta file.")
                .switch(
                    "description",
                    "parse the fasta header description",
                    Some('d'),
                )
                .category(Category::Experimental),
            Signature::build("from fa")
                .usage("Parse a fasta file.")
                .switch(
                    "description",
                    "parse the fasta header description",
                    Some('d'),
                )
                .category(Category::Experimental),
            Signature::build("from fa.gz")
                .usage("Parse a gzipped fasta file.")
                .switch(
                    "description",
                    "parse the fasta header description",
                    Some('d'),
                )
                .category(Category::Experimental),
            Signature::build("from fastq")
                .usage("Parse a fastq file.")
                .switch(
                    "description",
                    "parse the fastq header description",
                    Some('d'),
                )
                .switch(
                    "quality-scores",
                    "parse the fastq quality scores",
                    Some('q'),
                )
                .category(Category::Experimental),
            Signature::build("from fastq.gz")
                .usage("Parse a gzipped fastq file.")
                .switch(
                    "description",
                    "parse the fastq header description",
                    Some('d'),
                )
                .switch(
                    "quality-scores",
                    "parse the fastq quality scores",
                    Some('q'),
                )
                .category(Category::Experimental),
            Signature::build("from fq")
                .usage("Parse a fastq file.")
                .switch(
                    "description",
                    "parse the fastq header description",
                    Some('d'),
                )
                .switch(
                    "quality-scores",
                    "parse the fastq quality scores",
                    Some('q'),
                )
                .category(Category::Experimental),
            Signature::build("from fq.gz")
                .usage("Parse a gzipped fastq file.")
                .switch(
                    "description",
                    "parse the fastq header description",
                    Some('d'),
                )
                .switch(
                    "quality-scores",
                    "parse the fastq quality scores",
                    Some('q'),
                )
                .category(Category::Experimental),
            Signature::build("from bam")
                .usage("Parse a BAM file.")
                .category(Category::Experimental),
            Signature::build("from sam")
                .usage("Parse a SAM file.")
                .category(Category::Experimental),
            Signature::build("from cram")
                .usage("Parse a CRAM file into SAM output.")
                .category(Category::Experimental),
            Signature::build("from bcf")
                .usage("Parse a BCF file.")
                .category(Category::Experimental),
            Signature::build("from bcf.gz")
                .usage("Parse a gzipped BCF file.")
                .category(Category::Experimental),
            Signature::build("from vcf")
                .usage("Parse a VCF file.")
                .category(Category::Experimental),
            Signature::build("from vcf.gz")
                .usage("Parse a gzipped VCF file.")
                .category(Category::Experimental),
            Signature::build("from gff")
                .usage("Parse a GFF file.")
                .category(Category::Experimental),
        ]
    }

    fn run(
        &mut self,
        name: &str,
        call: &EvaluatedCall,
        input: &Value,
    ) -> Result<Value, LabeledError> {
        match name {
            "from fasta" => self.from_fasta(call, input, Compression::Uncompressed),
            "from fa" => self.from_fasta(call, input, Compression::Uncompressed),
            "from fastq" => self.from_fastq(call, input, Compression::Uncompressed),
            "from fq" => self.from_fastq(call, input, Compression::Uncompressed),
            "from fasta.gz" => self.from_fasta(call, input, Compression::Gzipped),
            "from fa.gz" => self.from_fasta(call, input, Compression::Gzipped),
            "from fastq.gz" => self.from_fastq(call, input, Compression::Gzipped),
            "from fq.gz" => self.from_fastq(call, input, Compression::Gzipped),
            "from bam" => self.from_bam(call, input),
            "from sam" => self.from_sam(call, input),
            "from cram" => self.from_cram(call, input),
            "from bcf" => self.from_bcf(call, input, Compression::Uncompressed),
            "from bcf.gz" => self.from_bcf(call, input, Compression::Gzipped),
            "from vcf" => self.from_vcf(call, input, Compression::Uncompressed),
            "from vcf.gz" => self.from_vcf(call, input, Compression::Gzipped),
            "from gff" => self.from_gff(call, input),
            _ => Err(LabeledError {
                label: "Plugin call with wrong name signature".into(),
                msg: "the signature used to call the plugin does not match any name in the plugin signature vector".into(),
                span: Some(call.head),
            }),
        }
    }
}
