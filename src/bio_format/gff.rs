/// The VCF format
use noodles::gff;
use nu_plugin::{EvaluatedCall, LabeledError};
use nu_protocol::Value;

use super::SpanExt;

/// The GFF3 headers
const GFF_COLUMNS: &[&str] = &[
    "ref_seq_name",
    "source",
    "ty",
    "start",
    "end",
    "score",
    "strand",
    "phase",
    "attributes",
];

/// Add a GFF record
fn add_record(call: &EvaluatedCall, r: gff::Record, vec_vals: &mut Vec<Value>) {
    let start = usize::from(r.start());
    let end = usize::from(r.end());

    let values_to_extend: Vec<Value> = vec![
        call.head.with_string(r.reference_sequence_name()),
        call.head.with_string(r.source()),
        call.head.with_string(r.ty()),
        Value::Int {
            val: start as i64,
            span: call.head,
        },
        Value::Int {
            val: end as i64,
            span: call.head,
        },
        call.head.with_string_or(r.score(), ""),
        call.head.with_string(r.strand()),
        call.head.with_string_or(r.phase(), ""),
        call.head.with_string(r.attributes()),
    ];

    vec_vals.extend_from_slice(&values_to_extend);
}

/// Parse a fasta file into a nushell structure.
pub fn from_gff_inner(call: &EvaluatedCall, input: &Value) -> Result<Vec<Value>, LabeledError> {
    // match on file type
    let stream = input.as_binary()?;

    let mut reader = gff::Reader::new(stream);

    let mut value_records = Vec::new();

    for record in reader.records() {
        let r = match record {
            Ok(rec) => rec,
            Err(e) => {
                return Err(LabeledError {
                    label: "Record reading failed.".into(),
                    msg: format!("cause of failure: {}", e),
                    span: Some(call.head),
                })
            }
        };

        let mut vec_vals = Vec::new();
        add_record(call, r, &mut vec_vals);

        value_records.push(Value::Record {
            cols: GFF_COLUMNS.iter().map(|e| String::from(*e)).collect(),
            vals: vec_vals,
            span: call.head,
        })
    }

    Ok(value_records)
}
