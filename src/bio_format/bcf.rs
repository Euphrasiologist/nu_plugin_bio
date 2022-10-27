/// The VCF format
use noodles::{
    bcf::{self, header::StringMaps},
    vcf,
};
use nu_plugin::{EvaluatedCall, LabeledError};
use nu_protocol::Value;

/// VCF column headers
const VCF_COLUMNS: &'static [&str] = &[
    "chrom",
    "pos",
    "rlen",
    "qual",
    "id",
    "ref",
    "alt",
    "filter",
    "info",
    "genotypes",
];

/// Add a VCF record to the vector.
/// TODO: make data more structured, so less is turned into a string immediately.
fn add_record(call: &EvaluatedCall, r: vcf::Record, vec_vals: &mut Vec<Value>) {
    let chrom = r.chromosome().to_string();
    let pos = usize::from(r.position());
    let rlen = r.reference_bases().len();
    let qual = match r.quality_score() {
        Some(q) => q.to_string(),
        None => "".into(),
    };
    let id = r.ids().to_string();
    let reference = r.reference_bases().to_string();
    let alt = r.alternate_bases().to_string();
    let filter = match r.filters() {
        Some(f) => f.to_string(),
        None => "".into(),
    };
    let info = r.info().to_string();
    let genotypes = r.genotypes().to_string();

    let values_to_extend: Vec<Value> = vec![
        Value::String {
            val: chrom,
            span: call.head,
        },
        Value::Int {
            val: pos as i64,
            span: call.head,
        },
        Value::Int {
            val: rlen as i64,
            span: call.head,
        },
        Value::String {
            val: qual,
            span: call.head,
        },
        Value::String {
            val: id,
            span: call.head,
        },
        Value::String {
            val: reference,
            span: call.head,
        },
        Value::String {
            val: alt,
            span: call.head,
        },
        Value::String {
            val: filter,
            span: call.head,
        },
        Value::String {
            val: info,
            span: call.head,
        },
        Value::String {
            val: genotypes,
            span: call.head,
        },
    ];

    vec_vals.extend_from_slice(&values_to_extend);
}

/// Parse a fasta file into a nushell structure.
pub fn from_bcf_inner(call: &EvaluatedCall, input: &Value) -> Result<Vec<Value>, LabeledError> {
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

    let mut reader = bcf::Reader::new(std::io::Cursor::new(stream));
    reader.read_file_format().unwrap();

    let raw_header = reader.read_header().unwrap();

    let header: vcf::Header = raw_header.parse().unwrap();
    let string_maps: StringMaps = raw_header.parse().unwrap();

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

        let v_r = r.try_into_vcf_record(&header, &string_maps).unwrap();

        let mut vec_vals = Vec::new();
        add_record(call, v_r, &mut vec_vals);

        value_records.push(Value::Record {
            cols: VCF_COLUMNS.iter().map(|e| String::from(*e)).collect(),
            vals: vec_vals,
            span: call.head,
        })
    }

    Ok(value_records)
}

/// Parse a fasta file into a nushell structure.
pub fn from_vcf_inner(call: &EvaluatedCall, input: &Value) -> Result<Vec<Value>, LabeledError> {
    // match on file type
    let stream = input.as_binary().unwrap();

    let mut reader = vcf::Reader::new(stream);
    let raw_header = reader.read_header().unwrap();

    let header: vcf::Header = match raw_header.parse() {
        Ok(h) => h,
        Err(_) => vcf::Header::default(),
    };

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
            cols: VCF_COLUMNS.iter().map(|e| String::from(*e)).collect(),
            vals: vec_vals,
            span: call.head,
        })
    }

    Ok(value_records)
}
