/// The fasta format
use noodles::fasta;
use nu_plugin::EvaluatedCall;
use nu_protocol::{Config, Value};

/// Parse a fasta file into a nushell structure.
pub fn from_fasta_inner(call: &EvaluatedCall, input: &Value) -> Vec<Value> {
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
        let r = record.unwrap();
        let id = r.name();
        let seq = r.sequence().as_ref();

        let mut vec_vals = Vec::new();

        vec_vals.push(Value::String {
            val: id.to_string(),
            span: call.head,
        });

        if description {
            let d_op = r.description();
            let d = match d_op {
                Some(e) => e,
                None => "",
            };

            vec_vals.push(Value::String {
                val: d.to_string(),
                span: call.head,
            });
        }

        vec_vals.push(Value::String {
            val: std::string::String::from_utf8(seq.to_owned()).unwrap(),
            span: call.head,
        });

        value_records.push(Value::Record {
            cols: cols.clone(),
            vals: vec_vals,
            span: call.head,
        })
    }

    value_records
}
