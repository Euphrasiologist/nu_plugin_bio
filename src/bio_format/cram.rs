/// The CRAM format
use noodles::cram;
use noodles::sam;
use nu_plugin::{EvaluatedCall, LabeledError};
use nu_protocol::Value;

use crate::bio_format::bam::{add_record, BAM_COLUMNS};

/// Parse a CRAM file into a nushell structure.
pub fn from_cram_inner(call: &EvaluatedCall, input: &Value) -> Result<Vec<Value>, LabeledError> {
    // match on file type
    let stream = match input {
        Value::Binary { val, span: _ } => val,
        other => {
            return Err(LabeledError {
                label: "Input should be binary.".into(),
                msg: format!("requires binary input, got {}", other.get_type()),
                span: Some(call.head),
            })
        }
    };

    let mut reader = cram::Reader::new(std::io::Cursor::new(stream));
    reader.read_file_definition().unwrap();
    let header: sam::Header = match reader.read_file_header() {
        Ok(s) => match s.parse() {
            Ok(s_p) => s_p,
            Err(e) => {
                return Err(LabeledError {
                    label: "Parsing SAM header from CRAM file failed.".into(),
                    msg: format!("cause of failure: {}", e),
                    span: Some(call.head),
                })
            }
        },
        Err(e) => {
            return Err(LabeledError {
                label: "CRAM file header reading failed.".into(),
                msg: format!("cause of failure: {}", e),
                span: Some(call.head),
            })
        }
    };

    let mut value_records = Vec::new();

    while let Some(container) = reader.read_data_container().unwrap() {
        for slice in container.slices() {
            // FIXME: if these unwraps in this section are converted to match LabeledErrors as
            // above, some of my test crams hang indefinitely in parsing(?)
            let records = slice.records(container.compression_header()).unwrap();

            for r in records {
                let r = r.try_into_alignment_record(&header).unwrap();
                let mut vec_vals = Vec::new();
                add_record(call, r, &mut vec_vals);

                value_records.push(Value::Record {
                    cols: BAM_COLUMNS.iter().map(|e| String::from(*e)).collect(),
                    vals: vec_vals,
                    span: call.head,
                })
            }
        }
    }

    Ok(value_records)
}
