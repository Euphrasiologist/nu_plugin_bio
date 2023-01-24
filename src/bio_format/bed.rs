use noodles::bed;
use nu_plugin::{EvaluatedCall, LabeledError};
use nu_protocol::Value;

/// BED reader type
const BED_COLUMN_NUMBER: u8 = 3;

/// Columns in a BAM/SAM file
pub const BED_COLUMNS: &[&str] = &[
    // Mandatory, name of chromosome
    "chrom",
    // Mandatory, start position
    "chromStart",
    // Mandatory, end position
    "chromEnd",
];

pub fn from_bed_inner(call: &EvaluatedCall, input: Value) -> Result<Vec<Value>, LabeledError> {
    let bytes = match input.as_binary() {
        Ok(b) => b,
        Err(e) => {
            return Err(LabeledError {
                label: "Value conversion to binary failed.".into(),
                msg: format!("cause of failure: {}", e),
                span: Some(call.head),
            })
        }
    };

    let mut reader = bed::Reader::new(bytes);

    let mut records = Vec::new();

    for result in reader.records::<BED_COLUMN_NUMBER>() {
        let record = result.map_err(|e| LabeledError {
            label: "Failed reading a record in the BED file".into(),
            msg: format!("{e}"),
            span: Some(call.head),
        })?;

        let mut row = Vec::new();

        row.push(Value::String {
            val: record.reference_sequence_name().into(),
            span: call.head,
        });
        let start: usize = record.start_position().into();
        row.push(Value::Int {
            val: start as i64,
            span: call.head,
        });
        let end: usize = record.end_position().into();
        row.push(Value::Int {
            val: end as i64,
            span: call.head,
        });

        records.push(Value::Record {
            cols: BED_COLUMNS
                .iter()
                .map(|e| e.to_string())
                .collect::<Vec<String>>(),
            vals: row,
            span: call.head,
        });
    }

    Ok(records)
}
