use bstr::io::*;
/// The GFA format.
/// We'll concern ourselves only with version 1.0.
use gfa::{
    gfa::Line::*,
    optfields::{OptField, OptFieldVal},
    parser::GFAParser,
};
use nu_plugin::{EvaluatedCall, LabeledError};
use nu_protocol::Value;
use std::io::BufReader;

/// We do a lot of string conversion in this module,
/// so make a string from utf8 function with nice error
/// handling.
fn string_from_utf8(
    inner: Vec<u8>,
    call: &EvaluatedCall,
    context: &str,
) -> Result<String, LabeledError> {
    match String::from_utf8(inner) {
        Ok(s) => Ok(s),
        Err(e) => Err(LabeledError {
            label: "Could convert bytes to string.".into(),
            msg: format!("{}: {}", context, e),
            span: Some(call.head),
        }),
    }
}

/// Parse a string representation of the option fields, until
/// we can come up with some better parsing.
fn parse_optfieldval(opt_field: OptField, call: &EvaluatedCall) -> Result<Value, LabeledError> {
    let tag = opt_field.tag;
    let val = opt_field.value;

    // TAG:TYPE:VALUE
    let tag_type_value =
        |tag: [u8; 2], typ: String, value: String, b: String| -> Result<Value, LabeledError> {
            let tag_string = string_from_utf8(tag.to_vec(), call, "tag is malformed")?;

            Ok(Value::String {
                val: format!("{tag_string}:{typ}:{b}{value}"),
                span: call.head,
            })
        };

    match val {
        // A (character)
        OptFieldVal::A(a) => tag_type_value(
            tag,
            String::from("A"),
            string_from_utf8(vec![a], call, "'A' value malformed")?,
            "".into(),
        ),
        // i (integer)
        OptFieldVal::Int(i) => tag_type_value(tag, String::from("i"), i.to_string(), "".into()),
        // f (real number)
        OptFieldVal::Float(f) => tag_type_value(tag, String::from("f"), f.to_string(), "".into()),
        // Z (string)
        OptFieldVal::Z(z) => tag_type_value(
            tag,
            String::from("Z"),
            string_from_utf8(z, call, "Z value malformed")?,
            "".into(),
        ),
        // what's J? should probably error out here.
        OptFieldVal::J(_j) => todo!(),
        // H (hexadecimal array)
        OptFieldVal::H(h) => tag_type_value(
            tag,
            String::from("Z"),
            h.iter()
                .map(|e| format!("{:#05x}", e))
                .fold(String::new(), |a, b| a + &b + ","),
            "".into(),
        ),
        // B (general array) - here it's split
        OptFieldVal::BInt(bi) => tag_type_value(
            tag,
            String::from("B"),
            bi.iter()
                .map(|e| e.to_string())
                .fold(String::new(), |a, b| a + &b + ","),
            "i:".into(),
        ),
        OptFieldVal::BFloat(bf) => tag_type_value(
            tag,
            String::from("B"),
            bf.iter()
                .map(|e| e.to_string())
                .fold(String::new(), |a, b| a + &b + ","),
            "f:".into(),
        ),
    }
}

pub fn from_gfa_inner(call: &EvaluatedCall, input: &Value) -> Result<Value, LabeledError> {
    let parser: GFAParser<Vec<u8>, Vec<OptField>> = GFAParser::new();

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

    let lines = BufReader::new(bytes).byte_lines();

    let mut header_nuon = Vec::new();
    let mut segments_nuon = Vec::new();
    let mut links_nuon = Vec::new();

    for line in lines {
        let line = match line {
            Ok(l) => l,
            Err(e) => {
                return Err(LabeledError {
                    label: "Could not read a line in the GFA.".into(),
                    msg: format!("cause of failure: {}", e),
                    span: Some(call.head),
                })
            }
        };
        // if this not added then
        if line.is_empty() {
            continue;
        }

        match parser.parse_gfa_line(line.as_ref()) {
            Ok(parsed) => {
                // what sort of line do we have?
                match parsed {
                    Header(h) => {
                        let version = h
                            .version
                            .and_then(|e| String::from_utf8(e).ok())
                            .unwrap_or_else(|| "No version specified.".into());

                        let opts: Result<Vec<Value>, _> = h
                            .optional
                            .iter()
                            .map(|e| parse_optfieldval(e.clone(), call))
                            .collect();

                        header_nuon.push(Value::Record {
                            cols: vec!["version".into(), "optional_fields".into()],
                            vals: vec![
                                Value::String {
                                    val: version,
                                    span: call.head,
                                },
                                Value::List {
                                    vals: opts?,
                                    span: call.head,
                                },
                            ],
                            span: call.head,
                        })
                    }
                    Segment(s) => {
                        // parse as string
                        let name = string_from_utf8(s.name, call, "segment name malformed");
                        let opts: Result<Vec<Value>, _> = s
                            .optional
                            .iter()
                            .map(|e| parse_optfieldval(e.clone(), call))
                            .collect();
                        // parse as string
                        let seq = string_from_utf8(s.sequence, call, "segment sequence malformed")?;

                        segments_nuon.push(Value::Record {
                            cols: vec!["name".into(), "sequence".into(), "optional_fields".into()],
                            vals: vec![
                                Value::String {
                                    val: name?,
                                    span: call.head,
                                },
                                Value::String {
                                    val: seq,
                                    span: call.head,
                                },
                                Value::List {
                                    vals: opts?,
                                    span: call.head,
                                },
                            ],
                            span: call.head,
                        })
                    }
                    Link(l) => {
                        let fo = l.from_orient.to_string();
                        let to = l.to_orient.to_string();
                        let fs = string_from_utf8(l.from_segment, call, "from segment malformed");
                        let ts = string_from_utf8(l.to_segment, call, "to segment malformed");
                        let overlap =
                            string_from_utf8(l.overlap, call, "overlap (CIGAR) malformed");
                        let opts: Result<Vec<Value>, _> = l
                            .optional
                            .iter()
                            .map(|e| parse_optfieldval(e.clone(), call))
                            .collect();

                        links_nuon.push(Value::Record {
                            cols: vec![
                                "from_orient".into(),
                                "to_orient".into(),
                                "from_segment".into(),
                                "to_segment".into(),
                                "overlaps".into(),
                                "optional_fields".into(),
                            ],
                            vals: vec![
                                Value::String {
                                    val: fo,
                                    span: call.head,
                                },
                                Value::String {
                                    val: to,
                                    span: call.head,
                                },
                                Value::String {
                                    val: fs?,
                                    span: call.head,
                                },
                                Value::String {
                                    val: ts?,
                                    span: call.head,
                                },
                                Value::String {
                                    val: overlap?,
                                    span: call.head,
                                },
                                Value::List {
                                    vals: opts?,
                                    span: call.head,
                                },
                            ],
                            span: call.head,
                        })
                    }
                    // TODO:
                    Containment(_c) => todo!(),
                    Path(_p) => todo!(),
                }
            }
            // I don't have access to the .tolerance field...
            // Err(err) if err.can_safely_continue(&parser.tolerance) => (),
            Err(e) => {
                return Err(LabeledError {
                    label: "Could not stream input as binary.".into(),
                    msg: format!("cause of failure: {}", e),
                    span: Some(call.head),
                })
            }
        };
    }

    Ok(Value::Record {
        cols: vec!["header".into(), "segments".into(), "links".into()],
        vals: vec![
            header_nuon
                .first()
                .unwrap_or(&Value::String {
                    val: "No header.".into(),
                    span: call.head,
                })
                .clone(),
            Value::List {
                vals: segments_nuon,
                span: call.head,
            },
            Value::List {
                vals: links_nuon,
                span: call.head,
            },
        ],
        span: call.head,
    })
}
