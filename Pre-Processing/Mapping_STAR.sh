#!/bin/bash

export PATH=~/bin/:$PATH

# Define paths and settings 

REF_FILE=~/data/genome/Reference.fa
REF_DIR=~/data/genome/star/
FASTQ_DIR=~/data/fastq/
GTF=~/data/genome/GRCh38.gtf
THREADS=8

# Build STAR index (only need to do this once) 

STAR --runThreadN "$THREADS" --runMode genomeGenerate --genomeDir "$REF_DIR"  --genomeFastaFiles "$REF_FILE" --sjdbGTFfile "$GTF" --sjdbOverhang 129


for R1_FILE in "$FASTQ_DIR"/*_R1_trimmed.fastq

do   
	R2_FILE="${R1_FILE/_R1/_R2}"   
	SAMPLE_NAME="$(basename "$R1_FILE" _R1_trimmed.fastq)"
	BAM_FILE="$SAMPLE_NAME.bam"    
	
    # Alignment with STAR
	nice STAR --runThreadN "$THREADS" --genomeDir "$REF_DIR" --readFilesIn "$R1_FILE" "$R2_FILE" --outSAMtype BAM SortedByCoordinate --outFileNamePrefix "$SAMPLE_NAME" 
	
    # Remove unmapped reads; 
    # Keep only properly paired reads; 
    # Minimum mapping quality of 20;
	nice samtools view -@ "$THREADS" -F 0x04 -f 0x2 -q 20 -b -o "$BAM_FILE" "${SAMPLE_NAME}Aligned.sortedByCoord.out.bam"
	
    # Sort files
	nice samtools sort -@ "$THREADS" -o "$BAM_FILE" -T "$SAMPLE_NAME.tmp" "$BAM_FILE" 
    # Index files  
	nice samtools index -@ "$THREADS" "$BAM_FILE"   
done