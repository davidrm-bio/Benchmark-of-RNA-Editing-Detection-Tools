#!/bin/bash
export PATH=/binf-isilon/rennie/gsn480/scratch/bin/:$PATH 

# Step 2: Mapping - STAR Aligner
Reference=/binf-isilon/rennie/gsn480/data/scratch/refGenome/human/grch38_star
Duplicates=/binf-isilon/rennie/gsn480/scratch/bin/picard-tools/MarkDuplicates.jar

# Create SAM file + SAM to BAM
find . -name "*.fastq.gz" | sort | paste - - | while read readA readB ; do nice STAR --runThreadN 8 --readFilesIn $readA $readB --genomeDir $Reference --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${readA/_R1_trimmed.fastq.gz/aligned} ; done

# Filter: Only mapped, properly pair and with an alignment quality of 20
find . -name '*.bam' | sort | paste - | while read file ; do samtools view -F 0x04 -f 0x2 -q 20 -b $file -@ 5 -o ${file/.bam/_filt.bam} ; done 

# Sort BAM File
find . -name "*filt.bam" | sort | paste - | while read file ; do samtools sort $file -o ${file/.bam/_sort.bam} -@ 8 ; done

# Remove Duplicates
find . -name "*sort.bam" | sort | paste - | while read file ; do java -jar $Duplicates INPUT=$file OUTPUT=${file/.bam/rmdup.bam} METRICS_FILE=${file/.bam/duplication.info} REMOVE_DUPLICATES=true ; done
