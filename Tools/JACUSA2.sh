#!/bin/bash

# Add local tools to path
export PATH=~/bin/:$PATH
JACUSA2=/binf-isilon/rennie/gsn480/scratch/bin/JACUSA_v2.0.2-RC.jar

# Define paths and settings
Reference=~/data/genome/Reference.fa

# Add MD Tag
find . -name "*.bam" | sort | paste - | while read file ; do samtools calmd $file $Reference > ${file/.bam/_MD.bam} ; done

# JACUSA2
java -jar $JACUSA2 call-2 -a D -p 5 -r Jacusa.out  WT_mock_clone1.bam,WT_mock_clone2.bam,WT_mock_clone3.bam ADAR1KO_mock_clone1.bam,ADAR1KO_mock_clone2.bam,ADAR1KO_mock_clone3.bam