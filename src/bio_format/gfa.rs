use bstr::io::*;
/// The GFA format.
/// We'll concern ourselves only with version 1.0.
use gfa::{
    gfa::Line::*,
    optfields::{OptField, OptFieldVal},
    parser::GFAParser,
};
use nu_plugin::{EvaluatedCall, LabeledError};
use nu_protocol::{record, Value};
use std::io::{BufRead, BufReader};

use super::{Compression, SpanExt};
use noodles::bgzf;

/// Compression status of a VCF reader.
enum GFAReader<'a> {
    Uncompressed(bstr::io::ByteLines<std::io::BufReader<&'a [u8]>>),
    Compressed(bstr::io::ByteLines<bgzf::Reader<std::io::BufReader<&'a [u8]>>>),
}

/// We do a lot of string conversion in this module,
/// so make a string from utf8 function with nice error
/// handling.
fn string_from_utf8(
    inner: Vec<u8>,
    call: &EvaluatedCall,
    context: &str,
) -> Result<String, LabeledError> {
    String::from_utf8(inner).map_err(|e| LabeledError {
        label: "Could convert bytes to string.".into(),
        msg: format!("{}: {}", context, e),
        span: Some(call.head),
    })
}

/// Parse a string representation of the option fields, until
/// we can come up with some better parsing.
fn parse_optfieldval(opt_field: OptField, call: &EvaluatedCall) -> Result<Value, LabeledError> {
    let tag = opt_field.tag;
    let val = opt_field.value;

    // TAG:TYPE:VALUE
    let tag_type_value = |typ: String, value: String, b: String| -> Result<Value, LabeledError> {
        let tag_string = string_from_utf8(tag.to_vec(), call, "tag is malformed")?;

        Ok(call
            .head
            .with_string(format!("{tag_string}:{typ}:{b}{value}")))
    };

    match val {
        // A (character)
        OptFieldVal::A(a) => tag_type_value(
            "A".into(),
            string_from_utf8(vec![a], call, "'A' value malformed")?,
            "".into(),
        ),
        // i (integer)
        OptFieldVal::Int(i) => tag_type_value("i".into(), i.to_string(), "".into()),
        // f (real number)
        OptFieldVal::Float(f) => tag_type_value(String::from("f"), f.to_string(), "".into()),
        // Z (string)
        OptFieldVal::Z(z) => tag_type_value(
            String::from("Z"),
            string_from_utf8(z, call, "Z value malformed")?,
            "".into(),
        ),
        // J is JSON
        // just handle this as a string
        OptFieldVal::J(j) => tag_type_value(
            String::from("J"),
            string_from_utf8(j, call, "J JSON value malformed")?,
            "".into(),
        ),
        // H (hexadecimal array)
        OptFieldVal::H(h) => tag_type_value(
            String::from("Z"),
            h.iter()
                .map(|e| format!("{:#05x}", e))
                .fold(String::new(), |a, b| a + &b + ","),
            "".into(),
        ),
        // B (general array) - here it's split
        OptFieldVal::BInt(bi) => tag_type_value(
            String::from("B"),
            bi.iter()
                .map(|e| e.to_string())
                .fold(String::new(), |a, b| a + &b + ","),
            "i:".into(),
        ),
        OptFieldVal::BFloat(bf) => tag_type_value(
            String::from("B"),
            bf.iter()
                .map(|e| e.to_string())
                .fold(String::new(), |a, b| a + &b + ","),
            "f:".into(),
        ),
    }
}

/// Convert GFA byte lines to nuon, given a compression status.
#[allow(clippy::too_many_arguments)]
fn lines_to_nuon<R: BufRead>(
    gfa_reader: ByteLines<R>,
    parser: GFAParser<Vec<u8>, Vec<OptField>>,
    header_nuon: &mut Vec<Value>,
    segments_nuon: &mut Vec<Value>,
    links_nuon: &mut Vec<Value>,
    containments_nuon: &mut Vec<Value>,
    paths_nuon: &mut Vec<Value>,
    call: &EvaluatedCall,
) -> Result<(), LabeledError> {
    for line in gfa_reader {
        let line = line.map_err(|e| LabeledError {
            label: "Could not read a line in the GFA.".into(),
            msg: format!("cause of failure: {}", e),
            span: Some(call.head),
        })?;
        // if this not added then
        if line.is_empty() {
            continue;
        }

        match parser.parse_gfa_line(line.as_ref()) {
            Ok(parsed) => {
                // what sort of line do we have?
                match parsed {
                    Header(h) => {
                        let version = h.version.and_then(|e| String::from_utf8(e).ok());

                        let opts: Result<Vec<Value>, _> = h
                            .optional
                            .iter()
                            .map(|e| parse_optfieldval(e.clone(), call))
                            .collect();

                        header_nuon.push(Value::record (
                            record! {"version" => call.head.with_string_or(version, "No version specified"),
                            "optional_fields" => Value::list(opts?, call.head)},  call.head
                        ))
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

                        segments_nuon.push(Value::record(
                            record! {
                            "name" => call.head.with_string(name?),
                            "sequence" =>  call.head.with_string(seq),
                            "optional_fields" => Value::list(opts?, call.head),
                            },
                            call.head,
                        ))
                    }
                    Link(l) => {
                        let fs = string_from_utf8(l.from_segment, call, "from segment malformed")?;
                        let ts = string_from_utf8(l.to_segment, call, "to segment malformed")?;
                        let overlap =
                            string_from_utf8(l.overlap, call, "overlap (CIGAR) malformed")?;
                        let opts: Result<Vec<Value>, _> = l
                            .optional
                            .iter()
                            .map(|e| parse_optfieldval(e.clone(), call))
                            .collect();

                        links_nuon.push(Value::record(
                            record! {
                                "from_orient" => call.head.with_string(l.from_orient),
                                "to_orient" => call.head.with_string(l.to_orient),
                                "from_segment" => call.head.with_string(fs),
                                "to_segment" => call.head.with_string(ts),
                                "overlaps" => call.head.with_string(overlap),
                                "optional_fields" => Value::list(opts?, call.head),
                            },
                            call.head,
                        ))
                    }
                    Containment(c) => {
                        let containment_name =
                            string_from_utf8(c.contained_name, call, "containment name malformed");
                        let container_name =
                            string_from_utf8(c.container_name, call, "container name malformed");
                        let overlap =
                            string_from_utf8(c.overlap, call, "overlap (CIGAR) malformed");
                        let position = c.pos;
                        let opts: Result<Vec<Value>, _> = c
                            .optional
                            .iter()
                            .map(|e| parse_optfieldval(e.clone(), call))
                            .collect();

                        containments_nuon.push(Value::record(
                            record! {
                                "containment_name" => call.head.with_string(containment_name?),
                                "containment_orient" => call.head.with_string(c.contained_orient),
                                "container_name" => call.head.with_string(container_name?),
                                "container_orient" => call.head.with_string(c.container_orient),
                                "overlap" => call.head.with_string(overlap?),
                                "position" => Value::int(position as i64, call.head),
                                "optional_fields" => Value::list(opts?, call.head),
                            },
                            call.head,
                        ))
                    }
                    Path(p) => {
                        let path_name = string_from_utf8(p.path_name, call, "malformed path name");
                        let segment_names = string_from_utf8(
                            p.segment_names,
                            call,
                            "segment names in path malformed",
                        )?;
                        let overlaps: Vec<Value> = p
                            .overlaps
                            .iter()
                            .map(|e| call.head.with_string_or(e.as_ref(), ""))
                            .collect();
                        let opts: Result<Vec<Value>, LabeledError> = p
                            .optional
                            .iter()
                            .map(|e| parse_optfieldval(e.clone(), call))
                            .collect();

                        paths_nuon.push(Value::record(
                            record! {
                                "path_name" => call.head.with_string(path_name?),
                                "segment_names" => call.head.with_string(segment_names),
                                "overlaps" => Value::list(overlaps, call.head),
                                "optional_fields" => Value::list(opts?, call.head),
                            },
                            call.head,
                        ))
                    }
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
    Ok(())
}

pub fn from_gfa_inner(
    call: &EvaluatedCall,
    input: &Value,
    gz: Compression,
) -> Result<Value, LabeledError> {
    let parser: GFAParser<Vec<u8>, Vec<OptField>> = GFAParser::new();

    let bytes = input.as_binary().map_err(|e| LabeledError {
        label: "Value conversion to binary failed.".into(),
        msg: format!("cause of failure: {}", e),
        span: Some(call.head),
    })?;

    let reader = BufReader::new(bytes);
    let lines = match gz {
        Compression::Uncompressed => GFAReader::Uncompressed(reader.byte_lines()),
        Compression::Gzipped => GFAReader::Compressed(bgzf::Reader::new(reader).byte_lines()),
    };

    let mut header_nuon = Vec::new();
    let mut segments_nuon = Vec::new();
    let mut links_nuon = Vec::new();
    let mut containments_nuon = Vec::new();
    let mut paths_nuon = Vec::new();

    match lines {
        GFAReader::Uncompressed(ur) => lines_to_nuon(
            ur,
            parser,
            &mut header_nuon,
            &mut segments_nuon,
            &mut links_nuon,
            &mut containments_nuon,
            &mut paths_nuon,
            call,
        )?,
        GFAReader::Compressed(cr) => lines_to_nuon(
            cr,
            parser,
            &mut header_nuon,
            &mut segments_nuon,
            &mut links_nuon,
            &mut containments_nuon,
            &mut paths_nuon,
            call,
        )?,
    };

    Ok(Value::record(
        record! {
            "header" => header_nuon.first().unwrap_or(&call.head.with_string("No header")).clone(),
            "segments" => Value::list(segments_nuon, call.head),
            "links" => Value::list(links_nuon, call.head),
            "containments" => Value::list(containments_nuon, call.head),
            "paths" => Value::list(paths_nuon, call.head)
        },
        call.head,
    ))
}
