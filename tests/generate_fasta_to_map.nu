## generate all the test data starting from a small genome file.

print -e "Starting pipeline to generate test data."
print -e ""

print -e "Please make sure MBG, samtools, bgzip, minimap2, and bcftools are installed/in your PATH."
print -e ""

# samtools, bgzip, minimap2, bcftools
let samtools_version = $"(samtools --version | lines | first | str replace 'samtools ' '')"
let bgzip_version = $"(bgzip --version | lines | first | str replace 'bgzip \(htslib\) ' '')"
let minimap2_version = $"(~/minimap2/minimap2 --version)"
let bcftools_version = $"(bcftools --version | lines | first | str replace 'bcftools ' '')"
let mbg_version = $"(do -i { ~/MBG/bin/MBG --version } | complete | get stderr | lines | get 1)"

print -e $"Software versions:\nSAMTOOLS:\t($samtools_version)"
print -e $"BGZIP\t($bgzip_version)"
print -e --no-newline $"MINIMAP2\t($minimap2_version)"
print -e --no-newline $"\nBCFTOOLS\t($bcftools_version)"
print -e $"\nMBG\t($mbg_version)"
print -e ""

# delete old files
rm -f fasta_to_map.fa*
rm -f drAilAlti1.fa.fai
rm -f map*

# the genome
# we're using the fasta parser from nu_plugin_bio here below.
let genome = (open ./drAilAlti1.fa | get sequence)
# genome length
let genome_length = ($genome
    | str length 
    | get 0)

# sequences we will split the random genome subsequences on
let splitters = ["AACCTAG", "TTCGACT", "TTCGACT", "AACCAAG"]
# DNA bases to join the split sequences with.
let bases = ["A", "T", "G", "CAATTAA"]

print -e "Generating 100 random subsequences."

# get 100 random subsequences
for s in 1..100 {
    let ri_1 = (random integer 0..$genome_length)
    let ri_2 = (random integer $ri_1..$genome_length)

    # get the subsequence
    let subseq = ($genome | str substring $ri_1..$ri_2)

    # random splitter to use
    let split_index = (random integer (0..(($splitters | length) - 1)))
    let split_string = ($splitters | get $split_index)

    #                      split each subseq                  and join with a random base
    echo $">sequence-($s)\n($subseq | split row $split_string | str join ($bases | get $split_index))\n"
        | save fasta_to_map.fa --append
}

print -e "BGZipping the simulated fasta file."
# make the zipped fasta to test
bgzip -c ./fasta_to_map.fa out> ./fasta_to_map.fa.gz
print -e ""
print -e "Using minimap2 to map reads to reference."
# # map synthesised subreads to genome
~/minimap2/minimap2 -a ./drAilAlti1.fa ./fasta_to_map.fa.gz -o map.sam
print -e ""
print -e "Using MBG to generate a GFA fake MBG assembly."
~/MBG/bin/MBG -i ./fasta_to_map.fa -o map.gfa -k 501
bgzip -k map.gfa

print -e ""
print -e "Converting output SAM to BAM."
# convert sam to bam
samtools view -T ./drAilAlti1.fa -bS map.sam -o map.bam
print -e "Converting BAM to CRAM."
# convert bam to cram
samtools view -T ./drAilAlti1.fa -C -o map.cram map.bam
print -e "Sorting and indexing BAM."
# index bam
samtools sort map.bam -o map_sorted.bam
samtools index map_sorted.bam
print -e ""
print -e "Calling variants."
# call variants
bcftools mpileup -f ./drAilAlti1.fa map_sorted.bam -Ob -o temp.bcf 
bcftools call -vmO z -o map.vcf.gz temp.bcf
rm temp.bcf
# and also gzip this
bcftools view map.vcf.gz -o map.bcf
bgzip -k map.bcf