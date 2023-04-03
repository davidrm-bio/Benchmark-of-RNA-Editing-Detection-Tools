#!/bin/bash

export PATH=~/bin/:$PATH

#conda activate myenv


# Define paths and settings 
Reference=~/data/genome/Reference.fa

# Run REDItools
find . -name "*.bam" | sort | paste - | while read file ; do  python ~/bin/reditools2.0-master/src/cineca/reditools.py -S -C -f $file -r $Reference -o ${file/.bam/.output}  ; done