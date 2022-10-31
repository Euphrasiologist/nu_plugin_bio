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

# see the current file formats currently supported below
# now you can just use open, and the file extension will be auto-detected.
open ./tests/test.fasta
    | get id
```

The backend is a <a href="https://github.com/zaeleus/noodles/">`noodles`</a> wrapper, an excellent, all-Rust bioinformatics I/O library.

## Aims

Aim to support the following:
- [x] BAM 1.6
- [x] BCF 2.2
- [ ] BED
- [x] CRAM 3.0
- [x] FASTA
- [x] FASTQ
- [x] GFF3
- [ ] GTF 2.2
- [x] SAM 1.6
- [x] VCF 4.3

Plus maybe some generic functions such as:
- [ ] indexing
- [ ] viewing
- [ ] getting headers
- [ ] counting lengths
- [ ] changing between file formats