# Nushell bio

A bioinformatics plugin for nushell. This plugin parses most common bioinformatics formats into structured data so you can use them with nushell more effectively.

# Quick setup

Go and get nushell, it's great. I'm assuming you have the rust toolchain installed. Then come back!

```nu
# clone this repo
git clone https://github.com/Euphrasiologist/nu_plugin_bio
# change into the repo directory
cd nu_plugin_bio
# build
# it's quite a long compile time...
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
  - [x] bcf.gz 
- [x] VCF 4.3
  - [x] vcf.gz
- [x] BED(3 only right now)
- [x] CRAM 3.0
- [x] FASTA
  - [x] fa.gz 
- [x] FASTQ
  - [x] fq.gz
- [x] GFF3
- [ ] GTF 2.2
- [x] SAM 1.6
- [x] GFA 1.0
  - [x] gfa.gz

Note that performance will not be optimal with the current state of `nu_plugin`, as we cannot access the engine state of nushell, and therefore need to load entire data structures into memory. Testing still needs to be done on large files.

## More?

If there's a bioinformatics format you want to add, let me know, or add a PR.