#!/bin/bash
export PATH=/binf-isilon/rennie/gsn480/scratch/bin/:$PATH 
export PATH=/binf-isilon/rennie/gsn480/scratch/bin/hisat2-2.2.1/:$PATH

# Step 2: Mapping - HISAT2 Aligner

Reference=/binf-isilon/rennie/gsn480/data/scratch/refGenome/human/grch38_hisat2/genome

# Create SAM file + SAM to BAM
find . -name "*.fastq.gz" | sort | paste - - | while read readA readB ; do nice hisat2 -x $Reference  -1 $readA -2 $readB --summary-file ${readA/_R1_trimmed.fastq.gz/summary.txt}  -p 8 | samtools view -Sb  > ${readA/_R1_trimmed.fastq.gz/.bam} ; done