# Using `nu_plugin_bio`

## The `tests` directory

The tests data contains loads of files to test the parsing capabilities of `nu_plugin_bio`. For small to medium files, parsing should be quick, and then we can use all the commands we can use with `nushell`.

For files which can be encountered gzipped in the wild, `nu_plugin_bio` provides parsers for those too, though you have to parse them explicitly, e.g:

```nu
# parse gfa
open map.gfa
# this will open a binary stream! Uh-oh
open map.gfa.gz
# we have to explicitly parse
open map.gfa.gz | from gfa.gz
```

For example, we might want to see the header of a BCF file.

```nu
# plain old header record
open map.bcf | get header
# explore around
open map.bcf | get header.contig
open map.bcf | get header.contig.drAilAlti1.length # etc
```

Or look at the table of links in a GFA file.

```nu
open map.gfa | get links
```

And see where there are segments that overlap with more than 100 bases:

```nu
open map.gfa 
    | get links
    | each {|e| let cigar = ($e.overlaps | str replace 'M' ''); if ($cigar | into int) > 100 { echo $e } }
```

## Some *big* data

Let's download some data. This is the latest human reference protein annotations, coding sequences, protein, and RNA fasta files.

```console
curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000001405.40/download?include_annotation_type=GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA&filename=GCF_000001405.40.zip" -H "Accept: application/zip"
```

This is what the files look like. We got some big files.

```
───┬──────────────────────┬──────┬───────────┬─────────────
 # │         name         │ type │   size    │  modified   
───┼──────────────────────┼──────┼───────────┼─────────────
 0 │ cds_from_genomic.fna │ file │ 333.6 MiB │ 6 hours ago 
 1 │ genomic.gff          │ file │   1.3 GiB │ 6 hours ago 
 2 │ protein.faa          │ file │  96.5 MiB │ 6 hours ago 
 3 │ rna.fna              │ file │ 684.0 MiB │ 6 hours ago 
───┴──────────────────────┴──────┴───────────┴─────────────
```

We will look quickly at each of these files in turn. Firstly `cds_from_genomic.fna`.

```nu
benchmark { let cds = (open cds_from_genomic.fna | from fasta) }
# 16sec 936ms 32µs 325ns
```

It's not insanely fast to load the file, but now it's in memory, so things will hopefully be faster. We can now for example, find total sequence length:

```nu
$cds | each {|e| $e.sequence | str length } | math sum
```

Or do some string searching:

```nu
# seven results for this particular string.
$cds 
    | where sequence =~ 'TTTGGAGGCTGCATCGCTCAAATCTTCTTCATCCACGTCGTTGGTGGTGTGGAGATGGTGCTGCTCATAGCCATGGCCTTTGACAGATA'
    | length
```

That GFF is huge! Let's parse it and see how long it takes.

```nu
benchmark { let gff = (open genomic.gff) }
# 7min 12sec 743ms 933µs 562ns
```

Pretty long. And going through it to do anything is quite slow too at the moment.

After a few minutes we get `2186877`.