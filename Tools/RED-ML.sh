#!/bin/bash

Reference=~/data/genome/Reference_hg19.fa
dbsnp=~data/genome/snp138.txt
simplerepeat=~data/genome/hg19.simpleRepeat.bed
aluBed=~data/genome/hg19.alu.bed


find . -name "*.bam" | sort | paste - | while read file ; do perl ~/bin/RED-ML/bin/red_ML.pl --rnabam $file --reference $Reference --dbsnp $dbsnp --simpleRepeat $simplerepeat --alu $aluBed --outdir ${file/.bam/_output} -p 0.5 ; done



# OPTIONS USED

# rnabam           Sorted BAM file  (input)
# --dbsnp          A file containing the SNPs of the dbSNP138 
# --simpleRepeat   A file containing genome-wide simple repeat annotation in BED format
# --alu            A file containing genome-wide Alu repeat annotation in BED format
# --outdir         The output  directory
# --p              Detection threshold cutoff [Between 0 and 1 default is 0.5]
