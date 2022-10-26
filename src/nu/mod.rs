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
            Signature::build("from fa")
                .usage("Parse a fasta file.")
                .switch(
                    "description",
                    "parse the fasta header description",
                    Some('d'),
                )
                .category(Category::Experimental),
            Signature::build("from bam")
                .usage("Parse a BAM file.")
                .category(Category::Experimental),
            Signature::build("from sam")
                .usage("Parse a SAM file.")
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
            "from fasta" => self.from_fasta(call, input),
            "from fa" => self.from_fasta(call, input),
            "from bam" => self.from_bam(call, input),
            "from sam" => self.from_sam(call, input),
            _ => Err(LabeledError {
                label: "Plugin call with wrong name signature".into(),
                msg: "the signature used to call the plugin does not match any name in the plugin signature vector".into(),
                span: Some(call.head),
            }),
        }
    }
}
