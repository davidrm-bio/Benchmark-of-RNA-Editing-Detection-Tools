# Summary
This repository contains the supplementary material for the review article '_Benchmarking RNA editing detection tools_'. You can find several folders which contain the scripts employed in the study.

* [Pre-processing](https://github.com/davidrm-bio/Benchmark-of-RNA-Editing-Detection-Tools/tree/main/Pre-Processing): includes the Bash scripts used for the trimming, mapping and processing of FASTQ files using three different mappers.
* [Tools](https://github.com/davidrm-bio/Benchmark-of-RNA-Editing-Detection-Tools/tree/main/Tools): includes the Bash scripts used for calling each tool. These scripts include comments explaining the different options.
* [Downstream](https://github.com/davidrm-bio/Benchmark-of-RNA-Editing-Detection-Tools/tree/main/Downstream): includes the Python and R scripts used in the downstream processing of the data for generating the plots. 

This repository also includes the supplementary figures and tables of the review article. The supplementary tables can be downloaded from [here](https://github.com/davidrm-bio/Benchmark-of-RNA-Editing-Detection-Tools/blob/main/Supplementary_Tables.xlsx). These tables and the supplementary figures are also shown below. 

#### Supplementary Figure S1
![](https://github.com/davidrm-bio/Benchmark-of-RNA-Editing-Detection-Tools/blob/main/Supplementary_Figure_S1.png)
**Supplementary Figure S1.** Statistics for BCFtools. For each condition replicate samples (n = 3) were processed individually keeping RNA editing sites (RES) absent in the HEK293T cell line. The average and the standard error of the mean are reported for each one of the measurements included: A) and B) average number (#) of total RES; C) and D) average percentage (%) of RES in REDIportal; E) and F) average percentage of RES in Alu regions.  Different threshold values for the minor allele frequency (MAF) were chosen as shown above.

#### Supplementary Figure S2
![](https://github.com/davidrm-bio/Benchmark-of-RNA-Editing-Detection-Tools/blob/main/Supplementary_Figure_S2.png)
**Supplementary Figure S2.** Statistics for RED-ML. For each condition, replicate samples (n = 3) were processed individually keeping RNA editing sites (RES) absent in the HEK293T cell line. The average and standard error of the mean are reported for each one of the measurements included: A) and B) average number (#) of total RES; C) and D) average percentage (%) of RES in REDIportal; and E) and F) average percentage of RES in Alu.   Different values for the detection threshold were chosen as shown above.

#### Supplementary Figure S3
![](https://github.com/davidrm-bio/Benchmark-of-RNA-Editing-Detection-Tools/blob/main/Supplementary_Figure_S3.png)
**Supplementary Figure S3.** Statistics for REDItools2. For each condition, replicate samples (n = 3) were processed individually keeping RNA editing sites (RES) absent in the HEK293T cell line. The average and standard error of the mean are reported for each one of the measurements included: A) and B) average number (#) of total RES; C) and D) average percentage (%) of RES in REDIportal; E) and F) average percentage of RES in Alu regions. Different threshold values for the number of supporting reads were chosen as shown above.

#### Supplementary Figure S4
![](https://github.com/davidrm-bio/Benchmark-of-RNA-Editing-Detection-Tools/blob/main/Supplementary_Figure_S4.png)
**Supplementary Figure S4.** Statistics for SPRINT. For each condition, replicate samples (n = 3) were processed individually keeping RNA editing sites (RES) absent in the HEK293T cell line. The average and standard error of the mean are reported for each one of the measurements included: A) and B) average number (#) of total RES; C) and D) average percentage (%) of RES in REDIportal; E) and F) average percentage of RES in Alu regions.  Different threshold values for the number of supporting reads were chosen as shown above.

#### Supplementary Figure S5
![](https://github.com/davidrm-bio/Benchmark-of-RNA-Editing-Detection-Tools/blob/main/Supplementary_Figure_S5.png)
**Supplementary Figure S5.** Statistics for JACUSA2. For each condition, replicate samples (n = 3) were merged into one sample excluding RNA editing sites (RES) not present in the three replicates and present in the HEK293T cell line. Different measurements are reported: A) and B) total number (#) of RES; C) and D) percentage (%) of RES in REDIportal; E) and F) percentage of RES in Alu regions. Different threshold values for the number of supporting reads were chosen.

#### Supplementary Table S1
**Table S1**. Main features of 10 tools. The ease of installation section is based on our experience installing each tool and its required software products and tools. It is subjective, yet it should give an overview of how much effort and time are required to install a particular tool. The date of the latest version section is based on the date provided on the GitHub or website of each tool.

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
| **Date of the latest version** 	| Jul-21 	| February 2016 (v0.3.1) 	| Nov-19 	| August 2022 (v1.0) 	| July 2021 (v2.0.2) 	| November 2018 (v0.1.8) 	| Jul-22 	| May 2016 (v1.1) 	| February 2018 (v1.0) 	| Mar-18 	|

#### Supplementary Table S2
**Table S2**.  Number of reads per sample. The average number of reads refers to the average amount of reads among three RNA-seq aligners in the BAM files used for the analysis after performing the quality control with SAMtools. Only mapped and properly pair reads were kept with an alignment quality of 20 for further analysis. The accession number from the Sequence Read Archive (SRA) is indicated for each sample.

| **Sample** 	| **SRA #** 	| **Total # reads in Fastq files (in million reads)** 	| **Average # mapped reads in GRCh37 (in million reads)** 	| **Average # mapped in GRCh38 (in million reads)** 	|
|:---:	|:---:	|---:	|---:	|---:	|
| **WT clone 1** 	| SRR5564274 	| 272,136 	| 52.81 	| 36.75 	|
| **WT clone 2** 	| SRR5564275 	| 405,346 	| 82.28 	| 54.35 	|
| **WT clone 3** 	| SRR5564276 	| 309,569 	| 66.84 	| 47 	|
| **ADAR1-KO clone 1** 	| SRR5564272 	| 347,676 	| 75.99 	| 55.12 	|
| **ADAR1-KO clone 2** 	| SRR5564273 	| 328,446 	| 60.6 	| 39.5 	|
| **ADAR1-KO clone 3** 	| SRR5564268 	| 509,534 	| 94.89 	| 72.79 	|

#### Supplementary Table S3
**Table S3**. Statistics for BCFtools. The number of RES were obtained individually for each sample by removing RES present in the HEK293T cell line. The average number of RES for each condition and replicate is reported for different minor allele frequency threshold values.

| **Aligner** 	| **Sample Condition** 	| **Minor Allele Frequency** 	| **Average # RES** 	| **RES SEM** 	| **Average % RES in REDIportal** 	| **Average % RES in Alu** 	|
|:---:	|:---:	|:---:	|:---:	|:---:	|:---:	|:---:	|
| **BWA** 	| WT 	| 0 	| 181554 	| 10536.61 	| 42.09 	| 47.83 	|
| **BWA** 	| WT 	| 0.1 	| 35152 	| 4104.4 	| 70.89 	| 67.81 	|
| **BWA** 	| ADAR1KO 	| 0 	| 103625 	| 7315.45 	| 7.46 	| 12.26 	|
| **BWA** 	| ADAR1KO 	| 0.1 	| 10870 	| 1308.79 	| 11.15 	| 11.63 	|
| **HISAT2** 	| WT 	| 0 	| 94162 	| 4392.25 	| 51.65 	| 59.63 	|
| **HISAT2** 	| WT 	| 0.1 	| 18549 	| 1985.09 	| 78.22 	| 75.8 	|
| **HISAT2** 	| ADAR1KO 	| 0 	| 41978 	| 2707.34 	| 9.5 	| 18.53 	|
| **HISAT2** 	| ADAR1KO 	| 0.1 	| 4036 	| 419.46 	| 17.02 	| 21.08 	|
| **STAR** 	| WT 	| 0 	| 17977 	| 1629.96 	| 57.93 	| 55.94 	|
| **STAR** 	| WT 	| 0.1 	| 7973 	| 835.43 	| 76.13 	| 71.86 	|
| **STAR** 	| ADAR1KO 	| 0 	| 8410 	| 949.12 	| 6.98 	| 11.07 	|
| **STAR** 	| ADAR1KO 	| 0.1 	| 2010 	| 213.62 	| 9.42 	| 11.77 	|

#### Supplementary Table S4
**Table S4**. Statistics for RED-ML. The number of RES were obtained individually for each sample by removing RES present in the HEK293T cell line. The average number of RES for each condition and replicate is reported for different values for the detection threshold.

| **Aligner** 	| **Sample Condition** 	| **Detection Threshold** 	| **Average # RES** 	| **RES SEM** 	| **Average % RES in REDIportal** 	| **Average % RES in Alu** 	|
|:---:	|:---:	|---:	|---:	|---:	|---:	|---:	|
| **BWA** 	| WT 	| 0.5 	| 17110 	| 1400.42 	| 64.61 	| 60.94 	|
| **BWA** 	| WT 	| 0.6 	| 14699 	| 1223.83 	| 70.16 	| 66.88 	|
| **BWA** 	| WT 	| 0.7 	| 12117 	| 1009.75 	| 75.36 	| 72.77 	|
| **BWA** 	| WT 	| 0.8 	| 8767 	| 759.47 	| 80.1 	| 78.46 	|
| **BWA** 	| WT 	| 0.9 	| 3949 	| 334.52 	| 82.08 	| 80.46 	|
| **BWA** 	| ADAR1KO 	| 0.5 	| 7228 	| 674.88 	| 10.88 	| 7.87 	|
| **BWA** 	| ADAR1KO 	| 0.6 	| 5316 	| 476.35 	| 13 	| 9.81 	|
| **BWA** 	| ADAR1KO 	| 0.7 	| 3683 	| 298.57 	| 15.72 	| 12.26 	|
| **BWA** 	| ADAR1KO 	| 0.8 	| 2206 	| 186.55 	| 18.74 	| 14.59 	|
| **BWA** 	| ADAR1KO 	| 0.9 	| 962 	| 97.17 	| 21.59 	| 13.8 	|
| **HISAT2** 	| WT 	| 0.5 	| 6158 	| 677.88 	| 88.2 	| 87.7 	|
| **HISAT2** 	| WT 	| 0.6 	| 6158 	| 677.88 	| 88.2 	| 87.7 	|
| **HISAT2** 	| WT 	| 0.7 	| 6158 	| 677.88 	| 88.2 	| 87.7 	|
| **HISAT2** 	| WT 	| 0.8 	| 5017 	| 544.06 	| 90.14 	| 90.28 	|
| **HISAT2** 	| WT 	| 0.9 	| 1976 	| 216.33 	| 92.94 	| 93.84 	|
| **HISAT2** 	| ADAR1KO 	| 0.5 	| 918 	| 68.72 	| 26.22 	| 31.32 	|
| **HISAT2** 	| ADAR1KO 	| 0.6 	| 918 	| 68.72 	| 26.22 	| 31.32 	|
| **HISAT2** 	| ADAR1KO 	| 0.7 	| 918 	| 68.72 	| 26.22 	| 31.32 	|
| **HISAT2** 	| ADAR1KO 	| 0.8 	| 652 	| 56.56 	| 29.51 	| 35.84 	|
| **HISAT2** 	| ADAR1KO 	| 0.9 	| 206 	| 28.62 	| 31.94 	| 38.66 	|
| **STAR** 	| WT 	| 0.5 	| 17309 	| 2003.06 	| 89.57 	| 91.66 	|
| **STAR** 	| WT 	| 0.6 	| 17309 	| 2003.06 	| 89.57 	| 91.66 	|
| **STAR** 	| WT 	| 0.7 	| 17309 	| 2003.06 	| 89.57 	| 91.66 	|
| **STAR** 	| WT 	| 0.8 	| 13386 	| 1544.89 	| 90.92 	| 93.29 	|
| **STAR** 	| WT 	| 0.9 	| 4262 	| 507.57 	| 93.03 	| 94.66 	|
| **STAR** 	| ADAR1KO 	| 0.5 	| 2267 	| 228.96 	| 36.78 	| 49.72 	|
| **STAR** 	| ADAR1KO 	| 0.6 	| 2267 	| 228.96 	| 36.78 	| 49.72 	|
| **STAR** 	| ADAR1KO 	| 0.7 	| 2267 	| 228.96 	| 36.78 	| 49.72 	|
| **STAR** 	| ADAR1KO 	| 0.8 	| 1571 	| 166.37 	| 39.05 	| 52.39 	|
| **STAR** 	| ADAR1KO 	| 0.9 	| 396 	| 51.28 	| 35.64 	| 47.15 	|

#### Supplementary Table S5
**Table S5**. Statistics for REDItools2. The number of RES were obtained individually for each sample by removing RES present in the HEK293T cell line. The average number of RES for each condition and replicate is reported for different numbers of supporting reads.

| **Aligner** 	| **Sample Condition** 	| **Number of Supporting Reads** 	| **Average # RES** 	| **RES SEM** 	| **Average % RES in REDIportal** 	| **Average % RES in Alu** 	|
|:---:	|:---:	|---:	|---:	|---:	|---:	|---:	|
| **BWA** 	| WT 	| 2 	| 344646 	| 27900.06 	| 32.98 	| 37.76 	|
| **BWA** 	| WT 	| 4 	| 201839 	| 18471.41 	| 30.99 	| 33.71 	|
| **BWA** 	| WT 	| 6 	| 124626 	| 11907.42 	| 30.37 	| 31.44 	|
| **BWA** 	| WT 	| 8 	| 76506 	| 7503.22 	| 31.11 	| 30.31 	|
| **BWA** 	| WT 	| 10 	| 42891 	| 4224.65 	| 34.7 	| 30.74 	|
| **BWA** 	| ADAR1KO 	| 2 	| 257445 	| 35146.1 	| 5.72 	| 11.78 	|
| **BWA** 	| ADAR1KO 	| 4 	| 161879 	| 25395.06 	| 5.83 	| 10.22 	|
| **BWA** 	| ADAR1KO 	| 6 	| 103153 	| 16782.15 	| 6.2 	| 9 	|
| **BWA** 	| ADAR1KO 	| 8 	| 63231 	| 9975.63 	| 6.95 	| 7.65 	|
| **BWA** 	| ADAR1KO 	| 10 	| 33648 	| 4353.35 	| 8.78 	| 5.51 	|
| **HISAT2** 	| WT 	| 2 	| 174481 	| 12403.45 	| 42.14 	| 46.68 	|
| **HISAT2** 	| WT 	| 4 	| 96400 	| 7906.89 	| 40.47 	| 43.29 	|
| **HISAT2** 	| WT 	| 6 	| 56114 	| 4882.07 	| 40.77 	| 42.6 	|
| **HISAT2** 	| WT 	| 8 	| 31656 	| 2966.5 	| 43.67 	| 44.61 	|
| **HISAT2** 	| WT 	| 10 	| 15083 	| 1607.95 	| 54.28 	| 53.43 	|
| **HISAT2** 	| ADAR1KO 	| 2 	| 110643 	| 14416.19 	| 7.39 	| 13.85 	|
| **HISAT2** 	| ADAR1KO 	| 4 	| 66136 	| 10483.14 	| 7.07 	| 12.35 	|
| **HISAT2** 	| ADAR1KO 	| 6 	| 39385 	| 6761.79 	| 7.07 	| 11.61 	|
| **HISAT2** 	| ADAR1KO 	| 8 	| 21492 	| 3762.11 	| 7.45 	| 11.3 	|
| **HISAT2** 	| ADAR1KO 	| 10 	| 8234 	| 1338.59 	| 9.5 	| 12.04 	|
| **STAR** 	| WT 	| 2 	| 246040 	| 18695.85 	| 41.63 	| 47.81 	|
| **STAR** 	| WT 	| 4 	| 135834 	| 11908.85 	| 40.74 	| 44.92 	|
| **STAR** 	| WT 	| 6 	| 79304 	| 7542.18 	| 41.44 	| 44.4 	|
| **STAR** 	| WT 	| 8 	| 45010 	| 4572.32 	| 44.63 	| 46.48 	|
| **STAR** 	| WT 	| 10 	| 21734 	| 2470.72 	| 55.47 	| 55.34 	|
| **STAR** 	| ADAR1KO 	| 2 	| 157254 	| 20891.59 	| 7.21 	| 16.15 	|
| **STAR** 	| ADAR1KO 	| 4 	| 92963 	| 14926.38 	| 7.02 	| 14.47 	|
| **STAR** 	| ADAR1KO 	| 6 	| 55287 	| 9656.7 	| 7.1 	| 13.58 	|
| **STAR** 	| ADAR1KO 	| 8 	| 30074 	| 5458.24 	| 7.42 	| 13.15 	|
| **STAR** 	| ADAR1KO 	| 10 	| 11651 	| 1959.67 	| 9.24 	| 13.8 	|


#### Supplementary Table S6
**Table S6**. Statistics for SPRINT. The number of RES were obtained individually for each sample by removing RES present in the HEK293T cell line. The average number of RES for each condition and replicate is reported for different numbers of supporting reads.
| **Aligner** 	| **Sample Condition** 	| **Number of Supporting Reads** 	| **Average # RES** 	| **RES SEM** 	| **Average % RES in REDIportal** 	| **Average % RES in Alu** 	|
|:---:	|:---:	|---:	|---:	|---:	|---:	|---:	|
| **BWA** 	| WT 	| 2 	| 27707 	| 4147.23 	| 96.19 	| 93.9 	|
| **BWA** 	| WT 	| 4 	| 6948 	| 1149.79 	| 98.08 	| 94.6 	|
| **BWA** 	| WT 	| 6 	| 3315 	| 585.57 	| 98.02 	| 95.28 	|
| **BWA** 	| WT 	| 8 	| 1981 	| 384.11 	| 98.09 	| 95.66 	|
| **BWA** 	| WT 	| 10 	| 1358 	| 272.5 	| 98.06 	| 96.21 	|
| **BWA** 	| ADAR1KO 	| 2 	| 919 	| 20.07 	| 84.93 	| 95.45 	|
| **BWA** 	| ADAR1KO 	| 4 	| 203 	| 11.86 	| 88.36 	| 96.39 	|
| **BWA** 	| ADAR1KO 	| 6 	| 106 	| 7.21 	| 88.17 	| 94.69 	|
| **BWA** 	| ADAR1KO 	| 8 	| 70 	| 3.18 	| 85.28 	| 92.65 	|
| **BWA** 	| ADAR1KO 	| 10 	| 50 	| 2.08 	| 82.81 	| 92.26 	|
| **HISAT2** 	| WT 	| 2 	| 6903 	| 1061.71 	| 96.95 	| 91.35 	|
| **HISAT2** 	| WT 	| 4 	| 2185 	| 364.41 	| 97.74 	| 92.69 	|
| **HISAT2** 	| WT 	| 6 	| 1078 	| 206.85 	| 97.9 	| 93.26 	|
| **HISAT2** 	| WT 	| 8 	| 630 	| 126.37 	| 97.75 	| 93.25 	|
| **HISAT2** 	| WT 	| 10 	| 415 	| 90.38 	| 97.64 	| 93.71 	|
| **HISAT2** 	| ADAR1KO 	| 2 	| 58 	| 8.69 	| 65.92 	| 65.47 	|
| **HISAT2** 	| ADAR1KO 	| 4 	| 16 	| 6.36 	| 91.11 	| 89.88 	|
| **HISAT2** 	| ADAR1KO 	| 6 	| 7 	| 2.89 	| 100 	| 100 	|
| **HISAT2** 	| ADAR1KO 	| 8 	| 3 	| 1.2 	| 100 	| 100 	|
| **HISAT2** 	| ADAR1KO 	| 10 	| 3 	| 0.88 	| 100 	| 100 	|

#### Supplementary Table S7
**Table S7**. Statistics for JACUSA2. The number of RES was obtained by merging replicates into one sample for the two conditions (i.e., WT and ADAR1-KO). Only RES present in the three replicates and not present in the HEK293T cell line are reported for different numbers of supporting reads.

| **Aligner** 	| **Sample Condition** 	| **Number of Supporting Reads** 	| **# RES** 	| **% RES in REDIportal** 	| **% RES in Alu** 	|
|:---:	|:---:	|---:	|---:	|---:	|---:	|
| **BWA** 	| WT 	| 2 	| 27388 	| 50.86 	| 47.84 	|
| **BWA** 	| WT 	| 4 	| 27388 	| 50.86 	| 47.84 	|
| **BWA** 	| WT 	| 6 	| 25365 	| 49.34 	| 46.5 	|
| **BWA** 	| WT 	| 8 	| 22238 	| 46.39 	| 43.7 	|
| **BWA** 	| WT 	| 10 	| 20048 	| 43.81 	| 41.28 	|
| **BWA** 	| ADAR1KO 	| 2 	| 18606 	| 26.32 	| 24.5 	|
| **BWA** 	| ADAR1KO 	| 4 	| 18606 	| 26.32 	| 24.5 	|
| **BWA** 	| ADAR1KO 	| 6 	| 17602 	| 25.47 	| 23.75 	|
| **BWA** 	| ADAR1KO 	| 8 	| 16067 	| 23.96 	| 22.18 	|
| **BWA** 	| ADAR1KO 	| 10 	| 14941 	| 22.47 	| 20.79 	|
| **HISAT2** 	| WT 	| 2 	| 13252 	| 67.97 	| 65.17 	|
| **HISAT2** 	| WT 	| 4 	| 13252 	| 67.97 	| 65.17 	|
| **HISAT2** 	| WT 	| 6 	| 12208 	| 66.68 	| 63.81 	|
| **HISAT2** 	| WT 	| 8 	| 10704 	| 64.13 	| 61.38 	|
| **HISAT2** 	| WT 	| 10 	| 9642 	| 61.85 	| 59.23 	|
| **HISAT2** 	| ADAR1KO 	| 2 	| 7739 	| 39.79 	| 38.75 	|
| **HISAT2** 	| ADAR1KO 	| 4 	| 7739 	| 39.79 	| 38.75 	|
| **HISAT2** 	| ADAR1KO 	| 6 	| 7341 	| 38.67 	| 37.53 	|
| **HISAT2** 	| ADAR1KO 	| 8 	| 6761 	| 36.49 	| 35.29 	|
| **HISAT2** 	| ADAR1KO 	| 10 	| 6342 	| 34.64 	| 33.44 	|
| **STAR** 	| WT 	| 2 	| 20455 	| 62.46 	| 60.29 	|
| **STAR** 	| WT 	| 4 	| 20455 	| 62.46 	| 60.29 	|
| **STAR** 	| WT 	| 6 	| 18860 	| 61 	| 58.93 	|
| **STAR** 	| WT 	| 8 	| 16380 	| 58.11 	| 56.2 	|
| **STAR** 	| WT 	| 10 	| 14681 	| 55.36 	| 53.5 	|
| **STAR** 	| ADAR1KO 	| 2 	| 12431 	| 34.7 	| 34.77 	|
| **STAR** 	| ADAR1KO 	| 4 	| 12431 	| 34.7 	| 34.77 	|
| **STAR** 	| ADAR1KO 	| 6 	| 11764 	| 33.65 	| 33.66 	|
| **STAR** 	| ADAR1KO 	| 8 	| 10751 	| 31.67 	| 31.65 	|
| **STAR** 	| ADAR1KO 	| 10 	| 10029 	| 29.59 	| 29.51 	|




