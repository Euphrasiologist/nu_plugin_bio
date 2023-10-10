use std::io::{BufRead, BufReader};

/// The fasta format
use noodles::{bgzf, fasta, fastq};
use nu_plugin::{EvaluatedCall, LabeledError};
use nu_protocol::Value;

use crate::bio_format::{Compression, SpanExt};

/// Compression status of a fastq reader.
enum FastqReader<'a> {
    Uncompressed(Box<fastq::Reader<&'a [u8]>>),
    Compressed(Box<fastq::Reader<BufReader<bgzf::Reader<&'a [u8]>>>>),
}

/// Compression status of a fasta reader.
enum FastaReader<'a> {
    Uncompressed(Box<fasta::Reader<&'a [u8]>>),
    Compressed(fasta::Reader<Box<bgzf::Reader<&'a [u8]>>>),
}

/// Iterate over the records of a reader that implements [`BufRead`].
fn iterate_fastq_records<R: BufRead>(
    mut reader: fastq::Reader<R>,
    call: &EvaluatedCall,
    value_records: &mut Vec<Value>,
    description: bool,
    quality_scores: bool,
    cols: Vec<String>,
) -> Result<(), LabeledError> {
    // iterate over the records.
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
        vec_vals.push(call.head.with_string_from_utf8(r.name()));

        if description {
            vec_vals.push(call.head.with_string_from_utf8(r.description()));
        }

        if quality_scores {
            vec_vals.push(call.head.with_string_from_utf8(r.quality_scores()));
        }

        vec_vals.push(call.head.with_string_from_utf8(r.sequence()));

        value_records.push(Value::Record {
            cols: cols.clone(),
            vals: vec_vals,
            span: call.head,
        })
    }

    Ok(())
}

pub fn from_fastq_inner(
    call: &EvaluatedCall,
    input: &Value,
    gz: Compression,
) -> Result<Vec<Value>, LabeledError> {
    // parse description flag.
    let description = call.has_flag("description");
    let quality_scores = call.has_flag("quality-scores");

    let bytes = input.as_binary().map_err(|e| LabeledError {
        label: "Value conversion to binary failed.".into(),
        msg: format!("cause of failure: {}", e),
        span: Some(call.head),
    })?;

    let reader = match gz {
        Compression::Uncompressed => FastqReader::Uncompressed(Box::new(fastq::Reader::new(bytes))),
        Compression::Gzipped => {
            let gz = bgzf::Reader::new(bytes);
            FastqReader::Compressed(Box::new(fastq::Reader::new(BufReader::new(gz))))
        }
    };

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

    match reader {
        FastqReader::Uncompressed(u) => iterate_fastq_records(
            *u,
            call,
            &mut value_records,
            description,
            quality_scores,
            cols,
        )?,
        FastqReader::Compressed(c) => iterate_fastq_records(
            *c,
            call,
            &mut value_records,
            description,
            quality_scores,
            cols,
        )?,
    };

    Ok(value_records)
}

fn iterate_fasta_records<R: BufRead>(
    mut reader: fasta::Reader<R>,
    call: &EvaluatedCall,
    value_records: &mut Vec<Value>,
    description: bool,
    cols: Vec<String>,
) -> Result<(), LabeledError> {
    // iterate over the records
    for record in reader.records() {
        let r = record.map_err(|e| LabeledError {
            label: "Record reading failed.".into(),
            msg: format!("cause of failure: {}", e),
            span: Some(call.head),
        })?;

        let mut vec_vals = Vec::new();

        vec_vals.push(call.head.with_string(r.name()));

        if description {
            vec_vals.push(call.head.with_string_or(r.description(), ""));
        }

        vec_vals.push(call.head.with_string_from_utf8(r.sequence().as_ref()));

        value_records.push(Value::Record {
            cols: cols.clone(),
            vals: vec_vals,
            span: call.head,
        })
    }
    Ok(())
}

/// Parse a fasta file into a nushell structure.
pub fn from_fasta_inner(
    call: &EvaluatedCall,
    input: &Value,
    gz: Compression,
) -> Result<Vec<Value>, LabeledError> {
    // parse description flag.
    let description = call.has_flag("description");

    let bytes = input.as_binary()?;

    let reader = match gz {
        Compression::Uncompressed => FastaReader::Uncompressed(Box::new(fasta::Reader::new(bytes))),
        Compression::Gzipped => {
            let gz = Box::new(bgzf::Reader::new(bytes));
            FastaReader::Compressed(fasta::Reader::new(gz))
        }
    };

    let cols = match description {
        false => vec!["id".to_string(), "sequence".to_string()],
        true => vec![
            "id".to_string(),
            "description".to_string(),
            "sequence".to_string(),
        ],
    };

    let mut value_records = Vec::new();

    match reader {
        FastaReader::Uncompressed(u) => {
            iterate_fasta_records(*u, call, &mut value_records, description, cols)?
        }
        FastaReader::Compressed(c) => {
            iterate_fasta_records(c, call, &mut value_records, description, cols)?
        }
    };

    Ok(value_records)
}
