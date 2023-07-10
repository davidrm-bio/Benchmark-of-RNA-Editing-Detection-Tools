#!/bin/bash

export PATH=~/bin/:$PATH

#conda activate myenv


# Define paths and settings 
Reference=~/data/genome/Reference.fa

# Run REDItools
find . -name "*.bam" | sort | paste - | while read file ; do  python ~/bin/reditools2.0-master/src/cineca/reditools.py -S -C -bq 20 -q 20 -f $file -r $Reference -o ${file/.bam/.output}; done


# OPTIONS USED

# -S           Only sites with edits will be included in the output
# -C           Strand correction. Once the strand has been inferred, only bases according to this strand will be selected.
# -bq 20       The minimum base quality. Base quality below this value are not included
# -q           The minimum read quality. Reads whose mapping quality is below this value will be discarded.
