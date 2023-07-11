# Summary
This repository includes the supplementary material included in the review article about the benchmarking analysis of RNA editing detection tools. Additionally, there are several folders which include scripts employed in the study.

## Pre-processing
This folder contains Bash scripts used for the trimming, mapping and processing of the FASTQ files.

## Tools
This folder contains  Bash scripts for calling each tool. It includes comments explaining the options used.

## Dowsntream
This folder contains Python and R scripts used in the downstream processing of the data and for generating the plots. 

## Supplementary Material

### Supplementary Figure S1
![](https://github.com/davidrm-bio/Benchmark-of-RNA-Editing-Detection-Tools/blob/main/Supplementary_Figure_1.png)
**Supplementary Figure S1.** Statistics for BCFtools. For each condition replicate samples (n = 3) were processed individually keeping RNA editing sites (RES) absent in the HEK293T cell line. The average and the standard error of the mean are reported for each one of the measurements included: A) and B) average number (#) of total RES; C) and D) average percentage (%) of RES in REDIportal; E) and F) average percentage of RES in Alu regions.  Different threshold values for the minor allele frequency (MAF) were chosen as shown above.

### Supplementary Figure S2
![](https://github.com/davidrm-bio/Benchmark-of-RNA-Editing-Detection-Tools/blob/main/Supplementary_Figure_2.png)
**Supplementary Figure S2.** Statistics for RED-ML. For each condition, replicate samples (n = 3) were processed individually keeping RNA editing sites (RES) absent in the HEK293T cell line. The average and standard error of the mean are reported for each one of the measurements included: A) and B) average number (#) of total RES; C) and D) average percentage (%) of RES in REDIportal; and E) and F) average percentage of RES in Alu.   Different values for the detection threshold were chosen as shown above.

### Supplementary Figure S3
![](https://github.com/davidrm-bio/Benchmark-of-RNA-Editing-Detection-Tools/blob/main/Supplementary_Figure_3.png)
**Supplementary Figure S3.** Statistics for REDItools2. For each condition, replicate samples (n = 3) were processed individually keeping RNA editing sites (RES) absent in the HEK293T cell line. The average and standard error of the mean are reported for each one of the measurements included: A) and B) average number (#) of total RES; C) and D) average percentage (%) of RES in REDIportal; E) and F) average percentage of RES in Alu regions. Different threshold values for the number of supporting reads were chosen as shown above.

### Supplementary Figure S4
![](https://github.com/davidrm-bio/Benchmark-of-RNA-Editing-Detection-Tools/blob/main/Supplementary_Figure_4.png)
**Supplementary Figure S4.** Statistics for SPRINT. For each condition, replicate samples (n = 3) were processed individually keeping RNA editing sites (RES) absent in the HEK293T cell line. The average and standard error of the mean are reported for each one of the measurements included: A) and B) average number (#) of total RES; C) and D) average percentage (%) of RES in REDIportal; E) and F) average percentage of RES in Alu regions.  Different threshold values for the number of supporting reads were chosen as shown above.

### Supplementary Figure S5
![](https://github.com/davidrm-bio/Benchmark-of-RNA-Editing-Detection-Tools/blob/main/Supplementary_Figure_5.png)
**Supplementary Figure S5.** Statistics for JACUSA2. For each condition, replicate samples (n = 3) were merged into one sample excluding RNA editing sites (RES) not present in the three replicates and present in the HEK293T cell line. Different measurements are reported: A) and B) total number (#) of RES; C) and D) percentage (%) of RES in REDIportal; E) and F) percentage of RES in Alu regions. Different threshold values for the number of supporting reads were chosen.

### Supplementary Table S1

| **Category** 	| **REDItools2** 	| **GIREMI** 	| **RES-Scanner2** 	| **RNAEditor** 	| **JACUSA2** 	| **SPRINT** 	| **RESIC** 	| **RDDpred** 	| **RED-ML** 	| **DeepRed** 	|
|:---:	|:---:	|:---:	|:---:	|:---:	|:---:	|:---:	|:---:	|:---:	|:---:	|:---:	|
| **Required software products and tools** 	| HTSlib, SAMtools, MPI implementation (In parallel version) 	| HTSlib, SAMtools, R 	| Python >= v3.3, Pysam, BWA, SOAPnuke, JAVA, PILON, BLAT 	| BWA, Picard tools (v1.119), GATK (v3.7), BLAT, Pysam, BEDtools, Python-Qt5, Matplotlib, NumPy, Java (v8) 	| JACUSA2Helper (R package) 	| BWA, SAMtools 	| Bowtie, SAMtools 	| SAMtools, BCFtools, BAMtools, WEKA 	| SAMtools, Perl modules 	| MATLAB 	|
| **Requested data format(s) as input file** 	| BAM 	| BAM + Filtered SNVs 	| FASTQ or BAM 	| FASTQ 	| BAM 	| FASTQ or BAM 	| FASTQ 	| BAM 	| BAM 	| List of SNVs 	|
| **Handling of strandness of RNA-seq reads** 	| Yes 	| Yes 	| Yes (dUTP protocol) 	| Yes, but no specified strand 	| Yes 	| Yes 	| Yes, but no specified strand 	| No 	| No 	| No 	|
| **Acceptancy of replicate data** 	| No 	| Yes 	| No 	| No 	| Yes 	| No 	| Yes 	| Yes 	| No 	| No 	|
| **Direct comparison of DNA- and RNA-seq data** 	| Yes 	| No 	| Yes 	| No 	| Yes 	| No 	| Yes 	| No 	| No 	| No 	|
| **Possibility of using only RNA-seq data** 	| Yes 	| Yes 	| No 	| Yes 	| Only with replicates 	| Yes 	| Yes 	| Yes 	| Yes 	| Yes 	|
| **De novo detection of RES** 	| Yes 	| Yes 	| Yes 	| Yes 	| Yes 	| Yes 	| Yes 	| Yes 	| Yes 	| Yes 	|
| **Inclusion of read mappers** 	| No 	| No 	| Optional (BWA) 	| Mandatory (BWA) 	| No 	| Optional (BWA) 	| Mandatory (Bowtie) 	| No 	| No 	| No 	|
| **Identification of single-nucleotide variant (SNV)** 	| Yes (Samtools) 	| No 	| Yes (Samtools) 	| Yes (GATK) 	| Yes 	| Yes 	| Yes 	| Yes 	| No 	| No 	|
| **Annotation of RES** 	| Yes 	| No 	| Yes 	| Yes 	| No 	| No 	| No 	| No 	| Yes (Within repeats) 	| No 	|
| **Filtering of the identified RES** 	| Filtering 	| P-value 	| Filtering + P-value 	| Filtering 	| Filtering 	| Filtering 	| Filtering 	| Likelihoods 	| Detection threshold 	| - 	|
| **Ease of installation** 	| Medium 	| Easy 	| Medium 	| Hard 	| Easy 	| Easy 	| Easy 	| Medium 	| Easy 	| Hard 	|
| **Date of the latest version** 	| Jul-21 	| February 2016 (v0.3.1) 	| Nov-19 	| August 2022 (v1.0) 	| July 2021 (v2.0.2) 	| November 2018 (v0.1.8) 	| Jul-22 	| May 2016 (v1.1) 	| Feb ruary 2018 (v1.0) 	| Mar-18 	|
----

The Supplementary Tables can be found [here](https://github.com/davidrm-bio/Benchmark-of-RNA-Editing-Detection-Tools/blob/main/Supplementary_Tables.xlsx).
