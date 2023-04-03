#!/bin/bash

# conda activate sprint_env

rmsk=~data/genome/hg38_repeat.bed
Reference=~data/genome/Reference.fa

find . -name "*.bam" | sort | paste - | while read file ; do ~/bin/SPRINT/bin/sprint_from_bam -rp $rmsk $file $Reference ${file/.bam/_output} ~/bin/samtools ; done