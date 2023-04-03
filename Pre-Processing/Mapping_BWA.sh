#!/bin/bash

export PATH=~/bin/:$PATH

# Define paths and settings

REF_FILE=~/data/genome/Reference.fa
FASTQ_DIR=~/data/fastq/
THREADS=10

bwa index "$REF_FiLE"

for R1_FILE in "$FASTQ_DIR"/*_R1_trimmed.fastq.gz

do
	R2_FILE="${R1_FILE/_R1/_R2}"
	SAMPLE_NAME="$(basename "$R1_FILE" _R1_trimmed.fastq.gz)"
	SAM_FILE="$SAMPLE_NAME.sam"
	BAM_FILE="$SAMPLE_NAME.bam"

    # Alignment with mem algorithm 
    
	nice bwa mem -t "$THREADS" "$REF_FILE" "$R1_FILE" "$R2_FILE" > "$SAM_FILE"

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