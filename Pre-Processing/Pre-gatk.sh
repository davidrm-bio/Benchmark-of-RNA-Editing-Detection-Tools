#!/bin/bash

export PATH=/binf-isilon/rennie/gsn480/scratch/bin/gatk-4.3.0.0/:$PATH

# Pre-processing
find . -name "*sortrmdup.bam" | sort | paste - | while read file ; do samtools addreplacerg -r '@RG\tID:samplename\tSM:samplename' -o ${file/.bam/.tmp.bam} $file -@ 5; done

find . -name "*tmp.bam" | sort | paste - | while read file ; do samtools index $file ; done
