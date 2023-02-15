#!/bin/bash

export PATH=/binf-isilon/rennie/gsn480/scratch/bin/:$PATH 

# Step 1: Trimming of fastq files. 
# Last 20 bases are trimmed. Output GZIP compress

find . -name "*.fastq" | sort | paste - | while read A; do fastx_trimmer -l 130 -Q33 -z -i $A -o ${A/.fastq/_trimmed.fastq.gz} ; done