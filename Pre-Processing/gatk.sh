#!/bin/bash

export PATH=/binf-isilon/rennie/gsn480/scratch/bin/gatk-4.3.0.0/:$PATH

# Uncomment the correct version:
# BWA: 
# Reference=/binf-isilon/rennie/gsn480/data/scratch/refGenome/human/grch38_bwa/GRCh38.primary_assembly.genome.fa.gz
# HISAT2: 
# Reference=/binf-isilon/rennie/gsn480/data/scratch/refGenome/human/grch38_hisat2/genome.fa
# STAR: 
# Reference=/binf-isilon/rennie/gsn480/data/scratch/refGenome/human/grch38.fasta

# Variant calling
find . -name "*sortrmdup.tmp.bam" | sort | paste - | while read file ; do gatk HaplotypeCaller -R $Reference -I $file -O ${file/_filt_sortrmdup.bam/_raw_variants.vcf} ; done

# Extract SNPs
find . -name "*variants.vcf" | sort | paste - | while read file ; do gatk SelectVariants -R $Reference -V $file -O ${file/_raw_variants.vcf/_raw_snps.vcf} -selectType SNP ; done

# Filter SNPs
find . -name "*raw_snps.vcf" | sort | paste - | while read file ; do gatk VariantFiltration -R $Reference -V $file -O ${file/_raw_snps.vcf/_filt_snps.vcf} -filter-name "QD_filter" -filter "QD < 2.0"  -filter-name "FS_filter" -filter "FS > 60.0" -filter-name "MQ_filter" -filter "MQ < 40.0" -filter-name "SOR_filter" -filter "SOR > 4.0" -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"; done

# Exclude filtered variants
find . -name "*_filt_snps.vcf" | sort | paste - | while read file ; do gatk SelectVariants -R $Reference -V $file -O ${file/_filt_snps.vcf/_PASS_snps.vcf} --exclude-filtered ; done