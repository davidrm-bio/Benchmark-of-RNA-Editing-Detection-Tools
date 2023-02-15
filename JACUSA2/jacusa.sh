#!/bin/bash

JACUSA2=/binf-isilon/rennie/gsn480/scratch/bin/JACUSA_v2.0.2-RC.jar

# WTm Vs AKOm
java -jar $JACUSA2 call-2 WT_mock_clone1_filt_sortrmdup.bam,WT_mock_clone2_filt_sortrmdup.bam,WT_mock_clone3_filt_sortrmdup.bam ADAR1KO_mock_clone1_filt_sortrmdup.bam,ADAR1KO_mock_clone2_filt_sortrmdup.bam,ADAR1KO_mock_clone3_filt_sortrmdup.bam  -r JACUSA_WTmAKOm.out -s

# WTm Vs AKO150m
java -jar $JACUSA2 call-2 WT_mock_clone1_filt_sortrmdup.bam,WT_mock_clone2_filt_sortrmdup.bam,WT_mock_clone3_filt_sortrmdup.bam ADAR1p150KO_mock_clone1_filt_sortrmdup.bam,ADAR1p150KO_mock_clone2_filt_sortrmdup.bam,ADAR1p150KO_mock_clone3_filt_sortrmdup.bam -r JACUSA_WTmAKO150m.out -s

# WTt Vs AKOt
java -jar $JACUSA2 call-2 WT_IFNb_clone1_filt_sortrmdup.bam,WT_IFNb_clone2_filt_sortrmdup.bam,WT_IFNb_clone3_filt_sortrmdup.bam ADAR1KO_IFNb_clone1_filt_sortrmdup.bam,ADAR1KO_IFNb_clone2_filt_sortrmdup.bam,ADAR1KO_IFNb_clone3_filt_sortrmdup.bam -r JACUSA_WTtAKOt.out -s

# WTt Vs AKO150t
java -jar $JACUSA2 call-2 WT_IFNb_clone1_filt_sortrmdup.bam,WT_IFNb_clone2_filt_sortrmdup.bam,WT_IFNb_clone3_filt_sortrmdup.bam ADAR1p150KO_IFNb_clone1_filt_sortrmdup.bam,ADAR1p150KO_IFNb_clone2_filt_sortrmdup.bam,ADAR1p150KO_IFNb_clone3_filt_sortrmdup.bam -r JACUSA_WTtAKO150t.out -s

# WTm Vs WTt
java -jar $JACUSA2 call-2 WT_mock_clone1_filt_sortrmdup.bam,WT_mock_clone2_filt_sortrmdup.bam,WT_mock_clone3_filt_sortrmdup.bam WT_IFNb_clone1_filt_sortrmdup.bam,WT_IFNb_clone2_filt_sortrmdup.bam,WT_IFNb_clone3_filt_sortrmdup.bam -r JACUSA_WTmWTt.out -s

# AKOm Vs AKOt
java -jar $JACUSA2 call-2 ADAR1KO_mock_clone1_filt_sortrmdup.bam,ADAR1KO_mock_clone2_filt_sortrmdup.bam,ADAR1KO_mock_clone3_filt_sortrmdup.bam ADAR1KO_IFNb_clone1_filt_sortrmdup.bam,ADAR1KO_IFNb_clone2_filt_sortrmdup.bam,ADAR1KO_IFNb_clone3_filt_sortrmdup.bam -r JACUSA_AKOmAKOt.out -s

# AKO150m Vs AKO150t
java -jar $JACUSA2 call-2 ADAR1p150KO_mock_clone1_filt_sortrmdup.bam,ADAR1p150KO_mock_clone2_filt_sortrmdup.bam,ADAR1p150KO_mock_clone3_filt_sortrmdup.bam ADAR1p150KO_IFNb_clone1_filt_sortrmdup.bam,ADAR1p150KO_IFNb_clone2_filt_sortrmdup.bam,ADAR1p150KO_IFNb_clone3_filt_sortrmdup.bam -r JACUSA_AKO150mAKO150t.out -s
