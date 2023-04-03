#!/bin/bash
export PATH=~/bin/:$PATH

# Define paths and settings
Reference=~/data/genome/Reference.fa

find . -name "*.bam" | sort | paste - | while read file ; do bcftools  mpileup -Ou --max-depth 10000 -q 20 -Q 20 -f $Reference $file | bcftools call -mv -O b -o ${file/.bam/.bcf} ; done

