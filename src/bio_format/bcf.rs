/// The VCF format
use noodles::{
    bcf::{self, header::StringMaps},
    vcf,
};
use nu_plugin::{EvaluatedCall, LabeledError};
use nu_protocol::Value;

/// VCF column headers
const VCF_COLUMNS: &[&str] = &[
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

/// VCF header columns
const HEADER_COLUMNS: &[&str] = &[
    "file_format",
    "info",
    "filter",
    "format",
    "alt_alleles",
    "assembly",
    "contig",
    "meta",
    "pedigree",
    "samples",
];

/// This parses the header of a V/BCF
fn parse_header(call: &EvaluatedCall, h: &vcf::Header) -> Value {
    let file_format = Value::String {
        val: h.file_format().to_string(),
        span: call.head,
    };
    let infos = h.infos();

    // add infos into a record structure
    let infos_nuon = Value::Record {
        cols: infos.keys().map(|e| e.to_string()).collect(),
        vals: infos
            .values()
            .map(|f| Value::Record {
                cols: vec!["number".into(), "type".into(), "description".into()],
                vals: vec![
                    Value::String {
                        val: f.number().to_string(),
                        span: call.head,
                    },
                    Value::String {
                        val: f.ty().to_string(),
                        span: call.head,
                    },
                    Value::String {
                        val: f.description().to_string(),
                        span: call.head,
                    },
                ],
                span: call.head,
            })
            .collect(),
        span: call.head,
    };

    // the filters
    let filters = h.filters();
    let filters_nuon = Value::Record {
        cols: filters.keys().map(|e| e.to_string()).collect(),
        vals: filters
            .values()
            .map(|f| Value::Record {
                cols: vec!["description".into()],
                vals: vec![Value::String {
                    val: f.description().to_string(),
                    span: call.head,
                }],
                span: call.head,
            })
            .collect(),
        span: call.head,
    };

    // the formats
    let formats = h.formats();
    let formats_nuon = Value::Record {
        cols: formats.keys().map(|e| e.to_string()).collect(),
        vals: formats
            .values()
            .map(|f| Value::Record {
                cols: vec!["number".into(), "type".into(), "description".into()],
                vals: vec![
                    Value::String {
                        val: f.number().to_string(),
                        span: call.head,
                    },
                    Value::String {
                        val: f.ty().to_string(),
                        span: call.head,
                    },
                    Value::String {
                        val: f.description().to_string(),
                        span: call.head,
                    },
                ],
                span: call.head,
            })
            .collect(),
        span: call.head,
    };

    // alternative alleles
    let alt_alleles = h.alternative_alleles();
    let alt_alleles_nuon = Value::Record {
        cols: alt_alleles.keys().map(|e| e.to_string()).collect(),
        vals: alt_alleles
            .values()
            .map(|f| Value::Record {
                cols: vec!["description".into()],
                vals: vec![Value::String {
                    val: f.description().to_string(),
                    span: call.head,
                }],
                span: call.head,
            })
            .collect(),
        span: call.head,
    };

    // assembly
    let assembly = Value::String {
        val: h
            .assembly()
            .unwrap_or("No reference assembly URL specified.")
            .to_string(),
        span: call.head,
    };

    // contigs
    let contigs = h.contigs();
    let contigs_nuon = Value::Record {
        cols: contigs.keys().map(|e| e.to_string()).collect(),
        vals: contigs
            .values()
            .map(|f| {
                let mut cols = vec!["length".into()];
                cols.extend(f.other_fields().keys().map(|e| e.to_string()));

                let mut vals = vec![Value::Int {
                    val: f.length().unwrap_or(0) as i64,
                    span: call.head,
                }];

                vals.extend(f.other_fields().values().map(|e| Value::String {
                    val: e.clone(),
                    span: call.head,
                }));

                Value::Record {
                    cols,
                    vals,
                    span: call.head,
                }
            })
            .collect(),
        span: call.head,
    };

    // metadata
    let meta = h.meta();
    let meta_nuon = Value::Record {
        cols: meta.keys().map(|e| e.to_string()).collect(),
        vals: meta
            .values()
            .map(|f| Value::Record {
                cols: f.values().to_vec(),
                vals: f
                    .values()
                    .iter()
                    .map(|e| Value::String {
                        val: e.clone(),
                        span: call.head,
                    })
                    .collect(),
                span: call.head,
            })
            .collect(),
        span: call.head,
    };

    // don't know how to parse Pedigrees currently?
    let pedigree_db = h.pedigree_db();
    let pedigree_nuon = Value::String {
        val: pedigree_db.unwrap_or("No pedigree database.").into(),
        span: call.head,
    };

    // sample names
    let sample_names = h.sample_names();
    let sample_names_nuon = Value::List {
        vals: sample_names
            .iter()
            .map(|e| Value::String {
                val: e.clone(),
                span: call.head,
            })
            .collect(),
        span: call.head,
    };

    // TODO: I've skipped other records for the moment.
    // return the big record
    Value::Record {
        cols: HEADER_COLUMNS.iter().map(|e| e.to_string()).collect(),
        vals: vec![
            file_format,
            infos_nuon,
            filters_nuon,
            formats_nuon,
            alt_alleles_nuon,
            assembly,
            contigs_nuon,
            meta_nuon,
            pedigree_nuon,
            sample_names_nuon,
        ],
        span: call.head,
    }
}

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
pub fn from_bcf_inner(call: &EvaluatedCall, input: &Value) -> Result<Value, LabeledError> {
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

    match reader.read_file_format() {
        Ok(_) => (),
        Err(e) => {
            return Err(LabeledError {
                label: "Could not read file format".into(),
                msg: format!("file format unreadable due to: {}", e),
                span: Some(call.head),
            })
        }
    };

    // TODO: remove this unwrap
    let raw_header = match reader.read_header() {
        Ok(e) => e,
        Err(e) => {
            return Err(LabeledError {
                label: "Could not read header.".into(),
                msg: format!("header unreadable due to {}", e),
                span: Some(call.head),
            })
        }
    };

    // TODO: remove this unwrap
    let header: vcf::Header = raw_header.parse().unwrap();
    let header_nuon = parse_header(call, &header);
    // TODO: remove this unwrap
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

        // TODO: remove this unwrap
        let v_r = r.try_into_vcf_record(&header, &string_maps).unwrap();

        let mut vec_vals = Vec::new();
        add_record(call, v_r, &mut vec_vals);

        value_records.push(Value::Record {
            cols: VCF_COLUMNS.iter().map(|e| String::from(*e)).collect(),
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

/// Parse a fasta file into a nushell structure.
pub fn from_vcf_inner(call: &EvaluatedCall, input: &Value) -> Result<Value, LabeledError> {
    // match on file type
    // TODO: remove this unwrap
    let stream = input.as_binary().unwrap();

    let mut reader = vcf::Reader::new(stream);
    // TODO: remove this unwrap
    let raw_header = reader.read_header().unwrap();

    let header: vcf::Header = match raw_header.parse() {
        Ok(h) => h,
        Err(_) => vcf::Header::default(),
    };
    let header_nuon = parse_header(call, &header);

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
