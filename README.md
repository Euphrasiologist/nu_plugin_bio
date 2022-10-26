# Nushell bio

A bioinformatics plugin for nushell. The aim initially is to create a bunch of parsers for all of the common bioinformatics file formats and take it from there!

# Quick setup

Go and get nushell, it's great. I'm assuming you have the rust toolchain installed. Then come back!

```bash
# clone this repo
git clone https://github.com/Euphrasiologist/nu_plugin_bio
# build
cargo build --release
# register the plugin
register nu_plugin_bio/target/release/nu_plugin_bio

# parse a fasta!
open ./tests/test.fasta
```

## Aims

Aim to support the following:
- [] BAM 1.6
- [] BCF 2.2
- [] BED
- [] BGZF
- [] CRAM 3.0
- [-] FASTA
- [] FASTQ
- [] GFF3
- [] GTF 2.2
- [] SAM 1.6
- [] tabix
- [] VCF 4.3

Plus maybe some generic functions such as:
- [] indexing
- [] viewing
- [] getting headers
- [] counting lengths
- [] changing between file formats