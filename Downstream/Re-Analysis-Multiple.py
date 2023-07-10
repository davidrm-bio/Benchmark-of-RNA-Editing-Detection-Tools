import os
import pandas as pd
import ujson as json
from MainFunctions import *

# Global variables
path = '..'
aligners = [['bwa', 'hisat2', 'star'],  # For the rest
            ['bwa', 'hisat2']] # For SPRINT
samples = ['clone1', 'clone2', 'clone3']
conditions = ['WT', 'ADAR1KO']
#categories = ['N_RES', 'N_REDIportal', 'N_Alu', 'SNPs', 'List_of_Res']
thresholds = [['2', '4', '6', '8', '10'],  # For the rest
             ['0.5','0.6', '0.7', '0.8', '0.9'], # For REDML
             ['0', '0.1']] # For BCFtools
categories = ['N_res', 'N_db', 'N_Alu', 'List_of_Res']
tools = ['BCFTools', 'RED-ML', 'SPRINT', 'REDItools2', 'JACUSA2']
dbpath = '/binf-isilon/rennie/gsn480/data/compare_tools/dataset1/dbRNA-Editing'

# Load Data
reditools = load_database(path, 'Data_REDItool2.json')
redml = load_database(path, 'Data_REDML.json')
sprint = load_database(path, 'Data_SPRINT.json')
bcftools = load_database(path, 'Data_BCFTools.json')

print ('Loading REDIportal & Alu GRCh38')
rediportal38 = load_database(dbpath,'REDIportal.json')
alu38 = load_database(dbpath,'Alu_GRCh38.json')
print ('Finished...')

# Re-Analysis
Generate_Json_Multiple(path, 'Data_REDItools2-Multiple.json', reditools,rediportal38, alu38, aligners[0], thresholds[0])
print ('REDItools2 DONE')
Generate_Json_Multiple(path, 'Data_SPRINT-Multiple.json', sprint, rediportal38, alu38, aligners[1], thresholds[0])
print('SPRINT DONE')
Generate_Json_Multiple(path, 'Data_BCFTools-Multiple.json', bcftools,rediportal38, alu38, aligners[0], thresholds[2])
print ('BCFTools DONE')
print('Realising  memory ...\nLoading REDIportal & Alu GRCh37')
rediportal37 = load_database(dbpath,'REDIportal_hg19.json')
alu37 = load_database(dbpath,'Alu_GRCh37.json')
print ('Finished')
Generate_Json_Multiple(path, 'Data_REDML-Multiple.json', redml,rediportal37, alu37, aligners[0], thresholds[1])
print('RED-ML DONE')




