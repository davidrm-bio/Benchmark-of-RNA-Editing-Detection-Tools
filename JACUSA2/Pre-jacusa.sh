#!/bin/bash

# Pre-Processing script - Index BAM file and add MD tag

# Add local tools to path
export PATH=/binf-isilon/rennie/gsn480/scratch/bin/:$PATH

# Uncomment the correct version:
# BWA: 
# Reference=/binf-isilon/rennie/gsn480/data/scratch/refGenome/human/grch38_bwa/GRCh38.primary_assembly.genome.fa.gz
# HISAT2: 
# Reference=/binf-isilon/rennie/gsn480/data/scratch/refGenome/human/grch38_hisat2/genome.fa
# STAR: 
# Reference=/binf-isilon/rennie/gsn480/data/scratch/refGenome/human/grch38.fasta

# Add MD Tag - Uncomment if MD Tag not present
# find . -name "*sortrmdup.bam" | sort | paste - | while read file ; do samtools calmd $file $Reference > ${file/.bam/.MD.bam} ; done

# Index BAM file
find . -name "*.bam" | sort | paste - | while read file ; do samtools index $file ; done