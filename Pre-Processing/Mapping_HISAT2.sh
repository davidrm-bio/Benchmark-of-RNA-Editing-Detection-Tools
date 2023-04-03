#!/bin/bash

export PATH=~/bin/:$PATH


# Define paths and settings 

REF_FILE=~/data/genome/Reference.fa
FASTQ_DIR=~/data/fastq/
SPLICE=~/data/Genome/genome.ss
EXON=~/data/Genome/genome.exon

THREADS=8


# Build HISAT2 index (only need to do this once) 

nice hisat2-build  --exon "$EXON" --ss "$SPLICE"  "$REF_FILE" "$REF_FILE"

for R1_FILE in "$FASTQ_DIR"/*_R1_trimmed.fastq

do   
	R2_FILE="${R1_FILE/_R1/_R2}"   
	SAMPLE_NAME="$(basename "$R1_FILE" _R1_trimmed.fastq)"
	SAM_FILE="$SAMPLE_NAME.sam"  
	BAM_FILE="$SAMPLE_NAME.bam" 

	# Alignment with hisat2
	nice hisat2 -p "$THREADS" -x "$REF_FILE" -1 "$R1_FILE" -2 "$R2_FILE" -S "$SAM_FILE"

    # Remove unmapped reads; 
    # Keep only properly paired reads; 
    # Minimum mapping quality of 20;
	nice samtools view -@ "$THREADS" -F 0x04 -f 0x2 -q 20 -b -o "$BAM_FILE" "$SAM_FILE"  

    # Sort files
	nice samtools sort -@ "$THREADS" -o "$BAM_FILE" -T "$SAMPLE_NAME.tmp" "$BAM_FILE" 

    # Index files
	nice samtools index -@ "$THREADS" "$BAM_FILE"   
	rm "$SAM_FILE" 
done