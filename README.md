# Nushell bio

A bioinformatics plugin for nushell. 

The aim initially is to create a bunch of parsers for all of the common bioinformatics file formats and take it from there.

# Quick setup

Go and get nushell, it's great. I'm assuming you have the rust toolchain installed. Then come back!

```nu
# clone this repo
git clone https://github.com/Euphrasiologist/nu_plugin_bio
# change into the repo directory
cd nu_plugin_bio
# build
cargo build --release
# register the plugin
register nu_plugin_bio/target/release/nu_plugin_bio

# see the current file formats currently supported below
# now you can just use open, and the file extension will be auto-detected.

# there are some test files in the tests/ dir.
open ./tests/test.fasta
    | get id

# if you want to add flags you have to explicitly use from <x>
# e.g. if you want descriptions in fasta files to be parsed.

open --raw ./tests/test.fasta 
    | from fasta -d
    | first
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
- [x] GFA 1.0

And their BGZIP counterparts where appropriate (.vcf.gz, .fasta.gz, etc).

Note that performance will not be optimal with the current state of `nu_plugin`, as we cannot access the engine state of nushell, and therefore need to load entire data structures into memory. Testing still needs to be done on large files.