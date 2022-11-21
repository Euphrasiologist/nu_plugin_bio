/// The fasta format
use noodles::{fasta, fastq};
use nu_plugin::{EvaluatedCall, LabeledError};
use nu_protocol::{Config, Value};

// TODO: add fastq here too.
pub fn from_fastq_inner(call: &EvaluatedCall, input: &Value) -> Result<Vec<Value>, LabeledError> {
    // parse description flag.
    let description = call.has_flag("description");
    let quality_scores = call.has_flag("quality-scores");

    let cat_string = input.into_string("", &Config::default());

    let mut reader = fastq::Reader::new(std::io::Cursor::new(cat_string));

    let cols = match (description, quality_scores) {
        (false, false) => vec!["id".to_string(), "sequence".to_string()],
        (true, false) => vec![
            "id".to_string(),
            "description".to_string(),
            "sequence".to_string(),
        ],
        (false, true) => vec![
            "id".to_string(),
            "quality_scores".to_string(),
            "sequence".to_string(),
        ],
        (true, true) => vec![
            "id".to_string(),
            "description".to_string(),
            "quality_scores".to_string(),
            "sequence".to_string(),
        ],
    };

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
        // TODO: remove this unwrap
        let id = std::str::from_utf8(r.name()).unwrap();
        let seq = r.sequence();

        let mut vec_vals = Vec::new();

        vec_vals.push(Value::String {
            val: id.to_string(),
            span: call.head,
        });

        if description {
            let d_op = r.description();
            // TODO: remove this unwrap
            let d = std::str::from_utf8(d_op).unwrap();

            vec_vals.push(Value::String {
                val: d.to_string(),
                span: call.head,
            });
        }

        if quality_scores {
            let q_op = r.quality_scores();
            // TODO: remove this unwrap
            let q = std::str::from_utf8(q_op).unwrap();

            vec_vals.push(Value::String {
                val: q.to_string(),
                span: call.head,
            });
        }

        vec_vals.push(Value::String {
            val: String::from_utf8(seq.to_owned()).unwrap(),
            span: call.head,
        });

        value_records.push(Value::Record {
            cols: cols.clone(),
            vals: vec_vals,
            span: call.head,
        })
    }

    Ok(value_records)
}

/// Parse a fasta file into a nushell structure.
pub fn from_fasta_inner(call: &EvaluatedCall, input: &Value) -> Result<Vec<Value>, LabeledError> {
    // parse description flag.
    let description = call.has_flag("description");

    let cat_string = input.into_string("", &Config::default());

    let mut reader = fasta::Reader::new(std::io::Cursor::new(cat_string));

    let cols = match description {
        false => vec!["id".to_string(), "sequence".to_string()],
        true => vec![
            "id".to_string(),
            "description".to_string(),
            "sequence".to_string(),
        ],
    };

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
        let id = r.name();
        let seq = r.sequence().as_ref();

        let mut vec_vals = Vec::new();

        vec_vals.push(Value::String {
            val: id.to_string(),
            span: call.head,
        });

        if description {
            let d_op = r.description();
            let d = d_op.unwrap_or("");

            vec_vals.push(Value::String {
                val: d.to_string(),
                span: call.head,
            });
        }

        vec_vals.push(Value::String {
            val: String::from_utf8(seq.to_owned()).unwrap(),
            span: call.head,
        });

        value_records.push(Value::Record {
            cols: cols.clone(),
            vals: vec_vals,
            span: call.head,
        })
    }

    Ok(value_records)
}
