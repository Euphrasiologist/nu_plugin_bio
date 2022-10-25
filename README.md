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

