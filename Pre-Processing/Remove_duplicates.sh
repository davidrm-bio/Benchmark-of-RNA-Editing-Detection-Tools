#!/bin/bash

Duplicates=/binf-isilon/rennie/gsn480/scratch/bin/picard-tools/MarkDuplicates.jar

find . -name "*.bam" | sort | paste - | while read file ; do java -jar $Duplicates INPUT=$file OUTPUT=${file/.bam/rmdup.bam} METRICS_FILE=${file/.bam/duplication.info} REMOVE_DUPLICATES=true ; done
