#!/bin/bash

# Export to PATH
export PATH=/binf-isilon/rennie/gsn480/scratch/bin/:$PATH
Duplicates=/binf-isilon/rennie/gsn480/scratch/bin/picard-tools/MarkDuplicates.jar

# Filter: Only mapped, properly pair and with an alignment quality of 20
find . -name '*.bam' | sort | paste - | while read file ; do samtools view -F 0x04 -f 0x2 -q 20 -b $file -@ 5 -o ${file/.bam/_filt.bam} ; done 

# Sort BAM File
find . -name "*filt.bam" | sort | paste - | while read file ; do samtools sort $file -o ${file/.bam/_sort.bam} -@ 8 ; done

# Remove Duplicates
find . -name "*sort.bam" | sort | paste - | while read file ; do java -jar $Duplicates INPUT=$file OUTPUT=${file/.bam/rmdup.bam} METRICS_FILE=${file/.bam/duplication.info} REMOVE_DUPLICATES=true ; done

