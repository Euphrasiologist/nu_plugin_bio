use nu_plugin::{EvaluatedCall, LabeledError};
use nu_protocol::Value;

use crate::bio_format::fasta::from_fasta_inner;
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
        let value_records = from_fasta_inner(call, input);
        Ok(Value::List {
            vals: value_records,
            span: call.head,
        })
    }
}
