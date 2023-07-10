#!/bin/bash
export PATH=~/bin/:$PATH

# Define paths and settings
Reference=~/data/genome/Reference.fa

find . -name "*.bam" | sort | paste - | while read file ; do bcftools  mpileup -Ou --max-depth 10000 -q 20 -Q 20 -f $Reference $file | bcftools call -mv -O b -o ${file/.bam/.bcf} ; done


# OPTIONS USED
# -Ou            output in uncompressed BCF format 
# --max-depth    Max raw per file depth
# -q 20          Skip alignments with a mapping quality score below 20
# -Q 20          Skip bases with a base quality score below 20
# -f $Reference  path to reference file

# -mv            Use model for multiallelic and rare-variant calling. Output only variants
# -O             Output in uncompressed BCF format
