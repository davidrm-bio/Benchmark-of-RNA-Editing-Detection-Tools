# Summary
This repository includes supplemental material for the benchmarking analysis of RNA editing detection tools. 

## Supplementary Figures & Tables Legends
The legend for the Supplementary Figures can be found below:

**Supplementary Figure S1.** Statistics for BCFtools. For each condition replicate samples (n = 3) were processed individually keeping RNA editing sites (RES) absent in the HEK293T cell line. The average and the standard error of the mean are reported for each one of the measurements included: A) and B) average number (#) of total RES; C) and D) average percentage (%) of RES in REDIportal; E) and F) average percentage of RES in Alu regions.  Different threshold values for the minor allele frequency (MAF) were chosen as shown above.

**Supplementary Figure S2.** Statistics for RED-ML. For each condition, replicate samples (n = 3) were processed individually keeping RNA editing sites (RES) absent in the HEK293T cell line. The average and standard error of the mean are reported for each one of the measurements included: A) and B) average number (#) of total RES; C) and D) average percentage (%) of RES in REDIportal; and E) and F) average percentage of RES in Alu.   Different values for the detection threshold were chosen as shown above.

**Supplementary Figure S3.** Statistics for REDItools2. For each condition, replicate samples (n = 3) were processed individually keeping RNA editing sites (RES) absent in the HEK293T cell line. The average and standard error of the mean are reported for each one of the measurements included: A) and B) average number (#) of total RES; C) and D) average percentage (%) of RES in REDIportal; E) and F) average percentage of RES in Alu regions. Different threshold values for the number of supporting reads were chosen as shown above.

**Supplementary Figure S4.** Statistics for SPRINT. For each condition, replicate samples (n = 3) were processed individually keeping RNA editing sites (RES) absent in the HEK293T cell line. The average and standard error of the mean are reported for each one of the measurements included: A) and B) average number (#) of total RES; C) and D) average percentage (%) of RES in REDIportal; E) and F) average percentage of RES in Alu regions.  Different threshold values for the number of supporting reads were chosen as shown above.

**Supplementary Figure S5.** Statistics for JACUSA2. For each condition, replicate samples (n = 3) were merged into one sample excluding RNA editing sites (RES) not present in the three replicates and present in the HEK293T cell line. Different measurements are reported: A) and B) total number (#) of RES; C) and D) percentage (%) of RES in REDIportal; E) and F) percentage of RES in Alu regions. Different threshold values for the number of supporting reads were chosen.

The legend for the Supplementary Tables can be found in the Excel file. Additionally, there are several folders that contain the scripts employed during the study.

## Pre-processing
This folder contains Bash scripts used for the trimming, mapping and processing of the FASTQ files.

## Tools
This folder contains  Bash scripts for calling each tool. It includes comments explaining the options used.

## Dowsntream
This folder contains Python and R scripts used in the downstream processing of the data and for generating the plots. 

