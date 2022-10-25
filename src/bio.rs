use std::vec;

use noodles::fasta;
use nu_plugin::{EvaluatedCall, LabeledError};
use nu_protocol::{Config, Value};
pub struct Bio;

impl Bio {
    // fn print_values(
    //     &self,
    //     _index: u32,
    //     _call: &EvaluatedCall,
    //     _input: &Value,
    // ) -> Result<(), LabeledError> {
    //     // Note. When debugging your plugin, you may want to print something to the console
    //     // Use the eprintln macro to print your messages. Trying to print to stdout will
    //     // cause a decoding error for your message
    //     // eprintln!("Calling test {} signature", index);
    //     // eprintln!("value received {:?}", input);

    //     Ok(())
    // }


    /// A minimal working example of parsing a fasta into Nushell.
    pub fn from_fasta(&self, call: &EvaluatedCall, input: &Value) -> Result<Value, LabeledError> {

        let cat_string = input.into_string("", &Config::default());

        let mut reader = fasta::Reader::new(std::io::Cursor::new(cat_string));

        let cols = vec!["id".to_string(), "sequence".to_string()];

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

        Ok(Value::List {
            vals: value_records,
            span: call.head,
        })
    }
}
