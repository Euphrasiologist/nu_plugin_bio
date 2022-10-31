use noodles::{
    bam,
    sam::{self, alignment::Record},
};
use nu_plugin::{EvaluatedCall, LabeledError};
use nu_protocol::Value;

/// Columns in a BAM/SAM file
pub const BAM_COLUMNS: &'static [&str] = &[
    "read_name",
    "flags",
    "reference_sequence_id",
    "alignment_start",
    "mapping_quality",
    "cigar",
    "mate_reference_sequence_id",
    "mate_alignment_start",
    "template_length",
    "sequence",
    "quality_scores",
    "data",
];

/// Parse a SAM record, and append to a vector
pub fn add_record(call: &EvaluatedCall, r: Record, vec_vals: &mut Vec<Value>) {
    let read_name = match r.read_name() {
        Some(r_n) => r_n.to_string(),
        None => "No read name.".into(),
    };

    let flags = r.flags().bits();
    let reference_sequence_id = match r.reference_sequence_id() {
        Some(r_s_id) => r_s_id.to_string(),
        None => "No reference sequence ID".into(),
    };
    let alignment_start = match r.alignment_start() {
        Some(a_s) => a_s.to_string(),
        None => "No alignment start".into(),
    };
    let mapping_quality = match r.mapping_quality() {
        Some(m_q) => format!("{}", u8::from(m_q)),
        None => "".to_string(),
    };
    let cigar = format!("{}", r.cigar());
    let mate_reference_sequence_id = match r.mate_reference_sequence_id() {
        Some(m_r_s) => format!("{}", m_r_s),
        None => "No mate reference sequence ID".into(),
    };
    let mate_alignment_start = match r.mate_alignment_start() {
        Some(m_a_s) => m_a_s.to_string(),
        None => "No mate alignment start".into(),
    };
    let template_length = r.template_length();
    let sequence: Vec<u8> = r.sequence().as_ref().iter().map(|e| u8::from(*e)).collect();
    let quality_scores = r.quality_scores().to_string();
    let data = r.data().to_string();

    let values_to_extend: Vec<Value> = vec![
        Value::String {
            val: read_name,
            span: call.head,
        },
        Value::String {
            val: format!("{:#06x}", flags),
            span: call.head,
        },
        Value::String {
            val: reference_sequence_id,
            span: call.head,
        },
        Value::String {
            val: alignment_start,
            span: call.head,
        },
        Value::String {
            val: mapping_quality,
            span: call.head,
        },
        Value::String {
            val: cigar,
            span: call.head,
        },
        Value::String {
            val: mate_reference_sequence_id,
            span: call.head,
        },
        Value::String {
            val: mate_alignment_start,
            span: call.head,
        },
        Value::Int {
            val: template_length.into(),
            span: call.head,
        },
        Value::String {
            val: std::string::String::from_utf8(sequence).unwrap(),
            span: call.head,
        },
        Value::String {
            val: quality_scores,
            span: call.head,
        },
        Value::String {
            val: data,
            span: call.head,
        },
    ];

    vec_vals.extend_from_slice(&values_to_extend);
}

/// Parse a BAM file into a nushell structure.
pub fn from_bam_inner(call: &EvaluatedCall, input: &Value) -> Result<Vec<Value>, LabeledError> {
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

    let mut reader = bam::Reader::new(std::io::Cursor::new(stream));
    let _ = reader.read_header();
    let _ = reader.read_reference_sequences();

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
            cols: BAM_COLUMNS.iter().map(|e| String::from(*e)).collect(),
            vals: vec_vals,
            span: call.head,
        })
    }

    Ok(value_records)
}

/// Parse a SAM file into a nushell structure.
pub fn from_sam_inner(call: &EvaluatedCall, input: &Value) -> Result<Vec<Value>, LabeledError> {
    // match on file type
    let stream = input.as_binary().unwrap();

    let mut reader = sam::Reader::new(stream);
    let header = reader.read_header().unwrap().parse().unwrap();

    let mut value_records = Vec::new();

    for record in reader.records(&header) {
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
            cols: BAM_COLUMNS.iter().map(|e| String::from(*e)).collect(),
            vals: vec_vals,
            span: call.head,
        })
    }

    Ok(value_records)
}
