#!/bin/bash

Reference="/binf-isilon/rennie/gsn480/data/compare_tools/dataset1/RED-ML/scratch/Genome/GRCh37.p13.genome.fa"
dbsnp="/binf-isilon/rennie/gsn480/data/compare_tools/dataset1/RED-ML/scratch/snp138.txt"
simplerepeat="/binf-isilon/rennie/gsn480/data/compare_tools/dataset1/RED-ML/scratch/hg19.simpleRepeat.bed"
aluBed="/binf-isilon/rennie/gsn480/data/compare_tools/dataset1/RED-ML/scratch/hg19.alu.bed"


find . -name "*.bam" | sort | paste - | while read file ; do perl /binf-isilon/rennie/gsn480/scratch/bin/RED-ML/bin/red_ML.pl --rnabam $file --reference $Reference --dbsnp $dbsnp --simpleRepeat $simplerepeat --alu $aluBed --outdir ${file/.bam/_output} -p 0.75 ; done