#!/bin/bash

# conda activate sprint_env

rmsk=~data/genome/hg38_repeat.bed
Reference=~data/genome/Reference.fa

find . -name "*.bam" | sort | paste - | while read file ; do ~/bin/SPRINT/bin/sprint_from_bam -rp $rmsk $file $Reference ${file/.bam/_output} ~/bin/samtools ; done

# OPTIONS USED

# -rp           Repear masker file obtained from USCS

# Default options used for cluster size.
# SPRINT cluster RES and only cluster with a count greater than the threshold
# are reported. For  Alu regions is 3/2; in non-Alu repeats is 5; outside
# repeat regions is 7. 
