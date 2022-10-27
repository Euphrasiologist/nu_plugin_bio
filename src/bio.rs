use crate::bio_format::bam::{from_bam_inner, from_sam_inner};
use crate::bio_format::bcf::{from_bcf_inner, from_vcf_inner};
use crate::bio_format::fasta::from_fasta_inner;
use nu_plugin::{EvaluatedCall, LabeledError};
use nu_protocol::Value;

pub struct Bio;

impl Bio {
    /// A minimal working example of parsing a fasta into Nushell.
    pub fn from_fasta(&self, call: &EvaluatedCall, input: &Value) -> Result<Value, LabeledError> {
        let value_records = from_fasta_inner(call, input)?;
        Ok(Value::List {
            vals: value_records,
            span: call.head,
        })
    }

    /// These B(S)AM functions are quite slow at the moment.
    pub fn from_bam(&self, call: &EvaluatedCall, input: &Value) -> Result<Value, LabeledError> {
        let value_records = from_bam_inner(call, input)?;
        Ok(Value::List {
            vals: value_records,
            span: call.head,
        })
    }
    /// These B(S)AM functions are quite slow at the moment.
    pub fn from_sam(&self, call: &EvaluatedCall, input: &Value) -> Result<Value, LabeledError> {
        let value_records = from_sam_inner(call, input)?;
        Ok(Value::List {
            vals: value_records,
            span: call.head,
        })
    }

    /// Parse a BCF.
    pub fn from_bcf(&self, call: &EvaluatedCall, input: &Value) -> Result<Value, LabeledError> {
        let value_records = from_bcf_inner(call, input)?;
        Ok(Value::List {
            vals: value_records,
            span: call.head,
        })
    }
    /// Parse a BCF.
    pub fn from_vcf(&self, call: &EvaluatedCall, input: &Value) -> Result<Value, LabeledError> {
        let value_records = from_vcf_inner(call, input)?;
        Ok(Value::List {
            vals: value_records,
            span: call.head,
        })
    }
}
