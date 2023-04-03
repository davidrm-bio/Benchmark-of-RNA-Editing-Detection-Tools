#!/bin/bash

export PATH=~/bin/:$PATH

# fastx_trimmer for trimming the last 20 nucleotides
find . -name "*.fastq" | sort | paste - | while read A; do fastx_trimmer -Q33 -l 130 -z -i $A -o ${A/.fastq/_trimmed.fastq.gz} ; done