use noodles::{
    bam,
    sam::{self, alignment::Record, header::record::value::Map},
};
use nu_plugin::{EvaluatedCall, LabeledError};
use nu_protocol::{Span, Value};

/// Columns in a BAM/SAM file
pub const BAM_COLUMNS: &[&str] = &[
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

/// Header fields in a B/SAM file
pub const HEADER_COLUMNS: &[&str] = &[
    "metadata",
    "reference_sequences",
    "read_groups",
    "programs",
    "comments",
];

trait SpanExt {
    fn string_value<S: ToString>(&self, s: S) -> Value;
    fn string_value_from_option<S: ToString>(&self, s: Option<S>, default: &str) -> Value;
}

impl SpanExt for Span {
    fn string_value<S: ToString>(&self, s: S) -> Value {
        Value::String {
            val: s.to_string(),
            span: *self,
        }
    }

    fn string_value_from_option<S: ToString>(&self, s: Option<S>, default: &str) -> Value {
        Value::String {
            val: s.map(|s| s.to_string()).unwrap_or(default.into()),
            span: *self,
        }
    }
}

/// Parse a B/SAM header
pub fn parse_header(call: &EvaluatedCall, h: &sam::Header) -> Value {
    // @HD in SAM.
    let header_op = h.header();
    // unwrap to the default header whatever that is?
    // should be able to modify it later to suit our needs.
    let default_map = Map::default();
    let header = header_op.unwrap_or(&default_map);

    let header_nuon = Value::Record {
        cols: vec![
            "version".into(),
            "sorting_order".into(),
            "grouping".into(),
            "sub_sort_order".into(),
        ],
        vals: vec![
            call.head.string_value(header.version()),
            // what's the default..?
            // if it's no good, we can always map -> string
            call.head
                .string_value(header.sort_order().unwrap_or_default()),
            // what's the default? see above.
            call.head
                .string_value(header.group_order().unwrap_or_default()),
            call.head
                .string_value_from_option(header.subsort_order(), "No subsort order."),
        ],
        span: call.head,
    };

    // @SQ.
    let reference_sequences = h.reference_sequences();
    let reference_sequences_nuon = Value::Record {
        cols: reference_sequences.keys().map(|e| e.to_string()).collect(),
        vals: reference_sequences
            .iter()
            .map(|(name, f)| Value::Record {
                cols: vec![
                    "sequence_name".into(),
                    "sequence_length".into(),
                    "alternate_locus".into(),
                    "alternate_names".into(),
                    "assembly_id".into(),
                    "description".into(),
                    "md5".into(),
                    "species".into(),
                    "molecule_topology".into(),
                    "uri".into(),
                ],
                vals: vec![
                    call.head.string_value(name),
                    Value::Int {
                        val: usize::from(f.length()) as i64,
                        span: call.head,
                    },
                    call.head
                        .string_value_from_option(f.alternative_locus(), "No alternative locus."),
                    call.head
                        .string_value_from_option(f.alternative_names(), "No alternative names."),
                    call.head
                        .string_value_from_option(f.assembly_id(), "No assembly ID."),
                    call.head
                        .string_value_from_option(f.description(), "No description"),
                    call.head
                        .string_value_from_option(f.md5_checksum(), "No md5 checksum."),
                    call.head
                        .string_value_from_option(f.species(), "No species name."),
                    call.head
                        .string_value_from_option(f.molecule_topology(), "No molecule topology."),
                    call.head.string_value_from_option(f.uri(), "No URI."),
                ],
                span: call.head,
            })
            .collect(),
        span: call.head,
    };

    // @RG
    let read_groups = h.read_groups();
    let read_groups_nuon = Value::Record {
        cols: read_groups.keys().cloned().collect(),
        vals: read_groups
            .iter()
            .map(|(id, f)| Value::Record {
                cols: vec![
                    "id".into(),
                    "barcode".into(),
                    "sequencing_center".into(),
                    "description".into(),
                    // no date?
                    // "date".into(),
                    "flow_order".into(),
                    "key_sequence".into(),
                    "library".into(),
                    "program".into(),
                    "predicted_insert_size".into(),
                    "platform".into(), // CAPILLARY, DNBSEQ (MGI/BGI), ELEMENT, HELICOS, ILLUMINA, IONTORRENT, LS454, ONT (Oxford Nanopore), PACBIO (Pacific Biosciences), SOLID, and ULTIMA
                    "platform_model".into(),
                    "platform_unit".into(),
                    "sample".into(),
                ],
                vals: vec![
                    call.head.string_value(id),
                    call.head
                        .string_value_from_option(f.barcode(), "No barcode."),
                    Value::String {
                        val: f
                            .sequencing_center()
                            .unwrap_or("No sequencing center.")
                            .into(),
                        span: call.head,
                    },
                    Value::String {
                        val: f.description().unwrap_or("No description.").into(),
                        span: call.head,
                    },
                    Value::String {
                        val: f.flow_order().unwrap_or("No flow order.").into(),
                        span: call.head,
                    },
                    Value::String {
                        val: f.key_sequence().unwrap_or("No key sequence.").into(),
                        span: call.head,
                    },
                    Value::String {
                        val: f.library().unwrap_or("No library.").into(),
                        span: call.head,
                    },
                    Value::String {
                        val: f.program().unwrap_or("No program.").into(),
                        span: call.head,
                    },
                    Value::Int {
                        val: f
                            .predicted_median_insert_size()
                            .map(|e| e as i64)
                            .unwrap_or(0),
                        span: call.head,
                    },
                    Value::String {
                        val: f
                            .platform()
                            .map(|e| e.to_string())
                            .unwrap_or_else(|| "No platform.".into()),
                        span: call.head,
                    },
                    call.head
                        .string_value_from_option(f.platform_model(), "No platform model"),
                    Value::String {
                        val: f.platform_unit().unwrap_or("No platform unit.").into(),
                        span: call.head,
                    },
                    Value::String {
                        val: f.sample().unwrap_or("No sample.").into(),
                        span: call.head,
                    },
                ],
                span: call.head,
            })
            .collect(),
        span: call.head,
    };

    // @PG
    let programs = h.programs();
    let programs_nuon = Value::Record {
        cols: programs.keys().cloned().collect(),
        vals: programs
            .iter()
            .map(|(id, f)| Value::Record {
                cols: vec![
                    "id".into(),
                    "name".into(),
                    "command_line".into(),
                    "previous_id".into(),
                    "description".into(),
                    "version".into(),
                ],
                vals: vec![
                    call.head.string_value(id),
                    call.head.string_value_from_option(f.name(), "No name."),
                    Value::String {
                        val: f
                            .command_line()
                            .map(|e| e.to_string())
                            .unwrap_or_else(|| "No command line.".into()),
                        span: call.head,
                    },
                    Value::String {
                        val: f
                            .previous_id()
                            .map(|e| e.to_string())
                            .unwrap_or_else(|| "No previous ID.".into()),
                        span: call.head,
                    },
                    Value::String {
                        val: f
                            .description()
                            .map(|e| e.to_string())
                            .unwrap_or_else(|| "No description.".into()),
                        span: call.head,
                    },
                    Value::String {
                        val: f
                            .version()
                            .map(|e| e.to_string())
                            .unwrap_or_else(|| "No version.".into()),
                        span: call.head,
                    },
                ],
                span: call.head,
            })
            .collect(),
        span: call.head,
    };

    // @CO
    let comments = h.comments();
    let comments_nuon = Value::List {
        vals: comments
            .iter()
            .map(|e| Value::String {
                val: e.clone(),
                span: call.head,
            })
            .collect(),
        span: call.head,
    };

    Value::Record {
        cols: HEADER_COLUMNS.iter().map(|e| e.to_string()).collect(),
        vals: vec![
            header_nuon,
            reference_sequences_nuon,
            read_groups_nuon,
            programs_nuon,
            comments_nuon,
        ],
        span: call.head,
    }
}

/// Parse a SAM record, and append to a vector
pub fn add_record(call: &EvaluatedCall, r: Record, vec_vals: &mut Vec<Value>) {
    let read_name = r
        .read_name()
        .map(|r_n| r_n.to_string())
        .unwrap_or("No read name.".into());

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
        call.head.string_value(read_name),
        call.head.string_value(format!("{:#06x}", flags)),
        call.head.string_value(reference_sequence_id),
        call.head.string_value(alignment_start),
        call.head.string_value(mapping_quality),
        call.head.string_value(cigar),
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
            val: String::from_utf8(sequence).unwrap(),
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
pub fn from_bam_inner(call: &EvaluatedCall, input: &Value) -> Result<Value, LabeledError> {
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

    let mut reader = bam::Reader::new(stream.as_slice());
    let raw_header = reader.read_header().map_err(|err| LabeledError {
        label: "Could not read header.".into(),
        msg: format!("error reading header at {}", err),
        span: Some(call.head),
    })?;

    // TODO: better error handling here.
    let header = if raw_header.is_empty() {
        let ref_seqs = reader.read_reference_sequences().unwrap();

        parse_header(
            call,
            &sam::Header::builder()
                .set_reference_sequences(ref_seqs)
                .build(),
        )
    } else {
        // this is required for reasons unclear to me...
        let _ = reader.read_reference_sequences().unwrap();
        parse_header(call, &raw_header.parse().unwrap())
    };

    let mut value_records = Vec::new();

    for record in reader.records() {
        let r = record.map_err(|e| LabeledError {
            label: "Record reading failed.".into(),
            msg: format!("cause of failure: {}", e),
            span: Some(call.head),
        })?;

        let mut vec_vals = Vec::new();

        add_record(call, r, &mut vec_vals);

        value_records.push(Value::Record {
            cols: BAM_COLUMNS.iter().map(|e| e.to_string()).collect(),
            vals: vec_vals,
            span: call.head,
        })
    }

    Ok(Value::Record {
        cols: vec!["header".into(), "body".into()],
        vals: vec![
            header,
            Value::List {
                vals: value_records,
                span: call.head,
            },
        ],
        span: call.head,
    })
}

/// Parse a SAM file into a nushell structure.
pub fn from_sam_inner(call: &EvaluatedCall, input: &Value) -> Result<Value, LabeledError> {
    // match on file type
    // TODO: remove this unwrap
    let stream = input.as_binary().unwrap();

    let mut reader = sam::Reader::new(stream);
    // TODO: remove this unwrap
    let header = reader.read_header().unwrap().parse().unwrap();
    let header_nuon = parse_header(call, &header);

    let mut value_records = Vec::new();

    for record in reader.records(&header) {
        let r = record.map_err(|e| LabeledError {
            label: "Record reading failed.".into(),
            msg: format!("cause of failure: {}", e),
            span: Some(call.head),
        })?;

        let mut vec_vals = Vec::new();
        add_record(call, r, &mut vec_vals);

        value_records.push(Value::Record {
            cols: BAM_COLUMNS.iter().map(|e| e.to_string()).collect(),
            vals: vec_vals,
            span: call.head,
        })
    }

    Ok(Value::Record {
        cols: vec!["header".into(), "body".into()],
        vals: vec![
            header_nuon,
            Value::List {
                vals: value_records,
                span: call.head,
            },
        ],
        span: call.head,
    })
}
