#!/bin/bash
export PATH=/binf-isilon/rennie/gsn480/scratch/bin/:$PATH 

# Step 2: Mapping - STAR Aligner

Reference=/binf-isilon/rennie/gsn480/data/scratch/refGenome/human/grch38_star

# Create SAM file + SAM to BAM
find . -name "*.fastq.gz" | sort | paste - - | while read readA readB ; do nice STAR --runThreadN 8 --readFilesIn $readA $readB --genomeDir $Reference --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${readA/_R1_trimmed.fastq.gz/aligned} ; done