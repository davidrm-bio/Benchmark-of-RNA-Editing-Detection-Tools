import os
import json
import pandas as pd
import warnings
import numpy as np
import time
from load_databases import generate_data

warnings.filterwarnings("ignore")

######################### Section 1: Load databases
start_time = time.time()

print ('\nLoading databases into memory, this will take time (~30 min)...')

dbSNP, REDIportal, Alu = generate_data()

print('Loading of dbSNP154, REDIportal & Alu finished. Time required {} min\n'.format (round ((time.time() - start_time)/60),2))

######################### Section 2: Analysis of the data

# Define variables
path = './'
files = os.listdir(path)
json_files = [file for file in files if 'json' in file]
wt_files = [file for file in json_files if 'WT' in file]
adar_files = [file for file in json_files if 'ADAR' in file]

# Load data
print ('Step 1: Load data')
wt_dict = dict()
for filename in wt_files:
    with open(filename) as json_file:
        file = json.load(json_file)
        for key in file:
            if key in wt_dict:
                wt_dict[key] +=1
            else:
                wt_dict[key] = 1
                
adar_dict = dict()
for filename in adar_files:
    with open(filename) as json_file:
        file = json.load(json_file)
        for key in file:
            if key in adar_dict:
                adar_dict[key] +=1
            else:
                adar_dict[key] = 1
                
# Retain only those present in at least 2 clones
## WT
wt_keys = list (wt_dict.keys())
wt_values = np.fromiter (wt_dict.values(),  dtype=float)
keep = list (np.where (wt_values >= 2)[0])
wt_keys_clean = [wt_keys[key] for key in keep]
wt_clean_tmp = dict(zip(wt_keys_clean, wt_keys_clean))

## ADAR
adar_keys = list (adar_dict.keys())
adar_values = np.fromiter (adar_dict.values(),  dtype=float)
keep = list (np.where (adar_values >= 2)[0])
adar_keys_clean = [adar_keys[key] for key in keep]
adar_clean_tmp = dict(zip(adar_keys_clean, adar_keys_clean))

# Remove variants from WT present in ADAR
wt_clean = [res for res in wt_clean_tmp if res not in adar_clean_tmp]
adar_clean = [res for res in adar_clean_tmp if res not in wt_clean_tmp]

# Remove possible SNPs
wt_clean_snp = [res for res in wt_clean if res not in dbSNP]
adar_clean_snp = [res for res in adar_clean if res not in dbSNP]

# Count how many A > G & C > T we have
A2G_wt = [ag for ag in wt_clean_snp if "AG" in ag]
T2C_wt = [ct for ct in wt_clean_snp if "TC" in ct]

A2G_adar = [ag for ag in adar_clean_snp if "AG" in ag]
T2C_adar = [ct for ct in adar_clean_snp if "TC" in ct]

A2G_T2C_wt = A2G_wt + T2C_wt
A2G_T2C_adar = A2G_adar + T2C_adar

# Check how many are in REDIportal
wt_db = [res for res in A2G_T2C_wt if res in REDIportal]
adar_db = [res for res in A2G_T2C_adar if res in REDIportal]

# FDR
G2A = [ga for ga in wt_clean_snp if "GA" in ga]
FDR = len(G2A) / (len(T2C_wt)+len(A2G_wt))

# Print out results

print ('Analysis with mpileup for 6 samples (x3 KO & x3 ADARKO)\n')
print ('There are {} SNVs in WT & {} in ADARKO\n'.format(len(wt_dict), len(adar_dict)))
print ('\nStep 2: Retain RES present in at least 2 clones')
print ('There are {} RES in WT & {} in ADAR1KO\n'.format(len(wt_clean_tmp), len(adar_clean_tmp)))
print ('\nStep 3: Remove WT RES present in ADAR1KO')
print ('There are {} RES in WT & {} in ADAR1KO\n'.format(len(wt_clean), len(adar_clean)))
print ('\nStep 4: Remove RES present in dbSNP154\n')
print ('There are {} RES in WT & {} in ADAR1KO'.format(len(wt_clean_snp), len(adar_clean_snp)))
print ('{} WT RES and {} ADAR1KO RES were removed\n'.format(len (wt_clean) - len (wt_clean_snp), len(adar_clean)-len(adar_clean_snp)))
print ('\nStep 5: Calculate % of RES present in REDIportal\n')
print ('{} % of WT RES and {} % of ADAR1KO RES are in REDIportal'.format(round (len(wt_db)/len(A2G_T2C_wt)*100,5), round (len(adar_db)/len(A2G_T2C_adar)*100,5)))
print ('\nStep 6: Calculate number of A > G & T > C RES\n')
print ('There are {} A > G & {} T > C in WT sample'.format(len (A2G_wt), len(T2C_wt)))
print ('\nStep 7: Calculate FDR G2A / (T2C+A2G)')
print('There are {} G>A and the FDR is {}'.format(len(G2A), FDR))

## Second part; check how many sites are in Alu repeats
cont = 0
for ID in A2G_T2C_wt:
    chrom = ID.split('|')[0]
    pos = ID.split('|')[1]
    if chrom in Alu:
        if pos in Alu[chrom]:
            cont += 1

print ('There are {} Percentage of IDs in Alu repeats'.format(round (cont/len(A2G_T2C_wt)*100,2)))

print ("Analysis finished")