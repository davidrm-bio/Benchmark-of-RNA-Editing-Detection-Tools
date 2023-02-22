#!/bin/bash

export PATH=/binf-isilon/rennie/gsn480/scratch/bin/:$PATH 
export PATH=/binf-isilon/rennie/gsn480/scratch/bin/hisat2-2.2.1/:$PATH

# Define paths and settings 

REF_FILE=/binf-isilon/rennie/gsn480/data/compare_tools/dataset1/RED-ML/scratch/Genome/hisat2/GRCh37.p13.genome.fa

FASTQ_DIR=/binf-isilon/rennie/gsn480/data/compare_tools/dataset1/RED-ML/hisat2/mapping

THREADS=8

# Build HISAT2 index (only need to do this once) 
nice hisat2-build "$REF_FILE" "$REF_FILE"


for R1_FILE in "$FASTQ_DIR"/*_R1_trimmed.fastq

do   
	R2_FILE="${R1_FILE/_R1/_R2}"   
	SAMPLE_NAME="$(basename "$R1_FILE" _R1_trimmed.fastq)"
	SAM_FILE="$SAMPLE_NAME.sam"  
	BAM_FILE="$SAMPLE_NAME.bam"    
	
	nice hisat2 -p "$THREADS" -x "$REF_FILE" -1 "$R1_FILE" -2 "$R2_FILE" -S "$SAM_FILE"
	nice samtools view -@ "$THREADS" -F 0x04 -f 0x2 -q 20 -b -o "$BAM_FILE" "$SAM_FILE"   
	nice samtools sort -@ "$THREADS" -o "$BAM_FILE" -T "$SAMPLE_NAME.tmp" "$BAM_FILE"   
	nice samtools index -@ "$THREADS" "$BAM_FILE"   
	rm "$SAM_FILE" 
done

# Mapping with Hisat2 to hg19
# Filter: Only mapped, properly pair and with an alignment quality of 20
# Sort
# Index