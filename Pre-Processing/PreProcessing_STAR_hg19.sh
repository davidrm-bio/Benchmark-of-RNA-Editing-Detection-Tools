#!/bin/bash

export PATH=/binf-isilon/rennie/gsn480/scratch/bin/:$PATH 

# Define paths and settings 

REF_FILE=/binf-isilon/rennie/gsn480/data/compare_tools/dataset1/RED-ML/scratch/Genome/GRCh37.p13.genome.fa

FASTQ_DIR=/binf-isilon/rennie/gsn480/data/compare_tools/dataset1/RED-ML/star/mapping

THREADS=8


for R1_FILE in "$FASTQ_DIR"/*_R1_trimmed.fastq

do   
	R2_FILE="${R1_FILE/_R1/_R2}"   
	SAMPLE_NAME="$(basename "$R1_FILE" _R1_trimmed.fastq)"
	BAM_FILE="$SAMPLE_NAME.bam"    
	
	nice STAR --runThreadN "$THREADS" --genomeDir /binf-isilon/rennie/gsn480/data/compare_tools/dataset1/RED-ML/scratch/Genome/star --readFilesIn "$R1_FILE" "$R2_FILE" \ --outSAMtype BAM SortedByCoordinate --outFileNamePrefix "$SAMPLE_NAME" 
	
	nice samtools view -@ "$THREADS" -F 0x04 -f 0x2 -q 20 -b -o "$BAM_FILE" "${SAMPLE_NAME}Aligned.sortedByCoord.out.bam"
	
	nice samtools sort -@ "$THREADS" -o "$BAM_FILE" -T "$SAMPLE_NAME.tmp" "$BAM_FILE"   
	nice samtools index -@ "$THREADS" "$BAM_FILE"   
done

# Mapping with STAR to hg19
# Filter: Only mapped, properly pair and with an alignment quality of 20
# Sort
# Index