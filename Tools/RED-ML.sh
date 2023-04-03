#!/bin/bash

Reference=~/data/genome/Reference_hg19.fa
dbsnp=~data/genome/snp138.txt
simplerepeat=~data/genome/hg19.simpleRepeat.bed
aluBed=~data/genome/hg19.alu.bed


find . -name "*.bam" | sort | paste - | while read file ; do perl ~/bin/RED-ML/bin/red_ML.pl --rnabam $file --reference $Reference --dbsnp $dbsnp --simpleRepeat $simplerepeat --alu $aluBed --outdir ${file/.bam/_output} -p 0.5 ; done