use crate::bio_format::SpanExt;
use noodles::{
    bam,
    sam::{self, alignment::Record as SAMRecord, header::record::value::Map},
};
use nu_plugin::{EvaluatedCall, LabeledError};
use nu_protocol::{record, Record, Value};

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

/// Parse a B/SAM header
pub fn parse_header(call: &EvaluatedCall, h: &sam::Header) -> Value {
    // @HD in SAM.
    let header_op = h.header();
    // unwrap to the default header whatever that is?
    // should be able to modify it later to suit our needs.
    let default_map = Map::default();
    let header = header_op.unwrap_or(&default_map);

    let header_nuon = Value::record(
        record!(
        "version" => call.head.with_string(header.version()),
        // what's the default..?
        // if it's no good, we can always map -> string
        "sorting_order" => call.head.with_string(header.sort_order().unwrap_or_default()),
        "grouping" => call.head.with_string(header.group_order().unwrap_or_default()),
        "sub_sort_order" => call.head.with_string_or(header.subsort_order(), "No subsort order.")
        ),
        call.head,
    );

    // @SQ.
    let reference_sequences = h.reference_sequences();
    let reference_sequences_nuon = Value::record(Record {
        cols: reference_sequences.keys().map(|e| e.to_string()).collect(),
        vals: reference_sequences
            .iter()
            .map(|(name, f)| Value::record(record! {
                    "sequence_name" => call.head.with_string(name),
                    "sequence_length" => Value::int(usize::from(f.length()) as i64, call.head),
                    "alternate_locus" => call.head.with_string_or(f.alternative_locus(), "No alternative locus."),
                    "alternate_names" => call.head.with_string_or(f.alternative_names(), "No alternative names."),
                    "assembly_id" => call.head.with_string_or(f.assembly_id(), "No assembly ID."),
                    "description" => call.head.with_string_or(f.description(), "No description"),
                    "md5" => call.head.with_string_or(f.md5_checksum(), "No md5 checksum"),
                    "species" => call.head.with_string_or(f.species(), "No species name"),
                    "molecule_topology" => call.head.with_string_or(f.molecule_topology(), "No molecule topology"),
                    "uri" => call.head.with_string_or(f.uri(), "No URI"),
                }, call.head)
            )
            .collect() },
        call.head
    );

    // @RG
    let read_groups = h.read_groups();
    let read_groups_nuon = Value::record(Record {
        cols: read_groups.keys().cloned().collect(),
        vals: read_groups
            .iter()
            .map(|(id, f)| Value::record(record! {
                    "id" => call.head.with_string(id),
                    "barcode" => call.head.with_string_or(f.barcode(), "No barcode"),
                    "sequencing_center" => call.head.with_string_or(f.sequencing_center(), "No sequencing center"),
                    "description" => call.head.with_string_or(f.description(), "No description"),
                    // no date?
                    // "date",
                    "flow_order" => call.head.with_string_or(f.flow_order(), "No flow order"),
                    "key_sequence" => call.head.with_string_or(f.key_sequence(), "No key sequence"),
                    "library" => call.head.with_string_or(f.library(), "No library"),
                    "program" => call.head.with_string_or(f.program(), "No program"),
                    "platform" => call.head.with_string_or(f.platform(), "No platform"),
                    "predicted_insert_size" => Value::int(
                         f.predicted_median_insert_size()
                                .map(|e| e as i64)
                                .unwrap_or(0),
                    call.head
                ), // CAPILLARY, DNBSEQ (MGI/BGI), ELEMENT, HELICOS, ILLUMINA, IONTORRENT, LS454, ONT (Oxford Nanopore), PACBIO (Pacific Biosciences), SOLID, and ULTIMA
                    "platform_model" => call.head.with_string_or(f.platform_model(), "No platform model"),
                    "platform_unit" => call.head.with_string_or(f.platform_unit(), "No platform unit"),
                    "sample" => call.head.with_string_or(f.sample(), "No sample"),
                }, call.head
            ))
            .collect() },
         call.head,
    );

    // @PG
    let programs = h.programs();
    let programs_nuon = Value::record(Record {
        cols: programs.keys().cloned().collect(),
        vals: programs
            .iter()
            .map(|(id, f)| Value::record(record! {
                    "id" => call.head.with_string(id),
                    "name" => call.head.with_string_or(f.name(), "No name"),
                    "command_line" => call.head.with_string_or(f.command_line(), "No command line"),
                    "previous_id" => call.head.with_string_or(f.previous_id(), "No previous ID"),
                    "description" => call.head.with_string_or(f.description(), "No description"),
                    "version" => call.head.with_string_or(f.version(), "No version"),
    }, call.head
            ))
            .collect() },
        call.head,
    );

    // @CO
    let comments = h.comments();
    let comments_nuon = Value::list(
        comments.iter().map(|e| call.head.with_string(e)).collect(),
        call.head,
    );

    Value::record(
        record! {
        HEADER_COLUMNS[0] => header_nuon,
        HEADER_COLUMNS[1] =>    reference_sequences_nuon,
        HEADER_COLUMNS[2] =>   read_groups_nuon,
        HEADER_COLUMNS[3] =>    programs_nuon,
        HEADER_COLUMNS[4] => comments_nuon,
        },
        call.head,
    )
}

/// Parse a SAM record, and append to a vector
pub fn create_record_values(call: &EvaluatedCall, r: SAMRecord) -> Vec<Value> {
    let flags = r.flags().bits();
    let mapping_quality = r
        .mapping_quality()
        .map(|m_q| format!("{}", u8::from(m_q)))
        .unwrap_or_default();
    let sequence: Vec<u8> = r.sequence().as_ref().iter().map(|e| u8::from(*e)).collect();

    vec![
        call.head.with_string_or(r.read_name(), "No read name."),
        call.head.with_string(format!("{:#06x}", flags)),
        call.head
            .with_string_or(r.reference_sequence_id(), "No reference sequence ID"),
        call.head
            .with_string_or(r.alignment_start(), "No alignment start"),
        call.head.with_string(mapping_quality),
        call.head.with_string(r.cigar()),
        call.head.with_string_or(
            r.mate_reference_sequence_id(),
            "No mate reference sequence ID",
        ),
        call.head
            .with_string_or(r.mate_alignment_start(), "No mate alignment start"),
        call.head.with_string(r.template_length()),
        call.head.with_string(String::from_utf8(sequence).unwrap()),
        call.head.with_string(r.quality_scores()),
        call.head.with_string(r.data()),
    ]
}

/// Parse a BAM file into a nushell structure.
pub fn from_bam_inner(call: &EvaluatedCall, input: &Value) -> Result<Value, LabeledError> {
    // match on file type
    let stream = match input {
        Value::Binary { val, .. } => val,
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

    let value_records = reader
        .records()
        .map(|record| {
            let r = record.map_err(|e| LabeledError {
                label: "Record reading failed.".into(),
                msg: format!("cause of failure: {}", e),
                span: Some(call.head),
            })?;

            let inner_record = Record::from_iter(
                BAM_COLUMNS
                    .iter()
                    .map(|e| e.to_string())
                    .zip(create_record_values(call, r)),
            );

            Ok(Value::record(inner_record, call.head))
        })
        .collect::<Result<Vec<_>, LabeledError>>()?;

    Ok(Value::record(
        record! {
            "header" => header,
            "body" => Value::list(value_records, call.head)
        },
        call.head,
    ))
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

    let value_records = reader
        .records(&header)
        .map(|record| {
            let r = record.map_err(|e| LabeledError {
                label: "Record reading failed.".into(),
                msg: format!("cause of failure: {}", e),
                span: Some(call.head),
            })?;

            let inner_record = Record::from_iter(
                BAM_COLUMNS
                    .iter()
                    .map(|e| e.to_string())
                    .zip(create_record_values(call, r)),
            );

            Ok(Value::record(inner_record, call.head))
        })
        .collect::<Result<Vec<_>, LabeledError>>()?;

    Ok(Value::record(
        record! {
            "header" => header_nuon,
            "body" => Value::list(value_records, call.head)
        },
        call.head,
    ))
}
