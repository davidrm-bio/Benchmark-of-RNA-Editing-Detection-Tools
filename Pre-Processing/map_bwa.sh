#!/bin/bash
export PATH=/binf-isilon/rennie/gsn480/scratch/bin/:$PATH 

# Step 2: Mapping - BWA Aligner

Reference=/binf-isilon/rennie/gsn480/data/scratch/refGenome/human/grch38_bwa/GRCh38.primary_assembly.genome.fa.gz

# bwa v0.7.7

# Create SAM file
find . -name "*.fastq.gz" | sort | paste - - | while read readA readB ; do nice bwa mem -t 8 -v3 $Reference $readA $readB > ${readA/_R1_trimmed.fastq.gz/.sam} ; done

# Convert SAM to BAM
find . -name "*.sam" | sort | paste - | while read file ; do samtools view -S -b $file > ${file/.sam/.bam}; done