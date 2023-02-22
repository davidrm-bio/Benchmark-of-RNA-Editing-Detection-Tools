#!/bin/bash

export PATH=/binf-isilon/rennie/gsn480/scratch/bin/:$PATH

# Activate correct conda env when running

# Define paths and settings 
Reference="/binf-isilon/rennie/gsn480/data/scratch/refGenome/human/grch38_bwa/GRCh38.primary_assembly.genome.fa.gz" 

# Run REDItools
find . -name "*.bam" | sort | paste - | while read file ; do  python /binf-isilon/rennie/gsn480/scratch/bin/reditools2.0-master/src/cineca/reditools.py -f $file -r $Reference -o ${file/.bam/.output} -C -t ${file/.bam/_tmp}  ; done