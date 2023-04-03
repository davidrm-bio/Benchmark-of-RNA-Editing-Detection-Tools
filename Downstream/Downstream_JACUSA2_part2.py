# Modules
import os
import pandas as pd
import numpy as np
import ujson as json
import warnings
import sys
import time
from load_databases import generate_data

warnings.filterwarnings("ignore")

# Functions
def load_file (filename):
    '''Load tab-delimited files into a DataFrame'''
    df = pd.read_csv(os.path.join(filename), delimiter="\t") 
    return df

######################### Section 1: Load databases
start_time = time.time()

print ('\nLoading databases into memory, this will take time (~30 min)...')

dbSNP, REDIportal, Alu = generate_data()

print('Loading of dbSNP154, REDIportal & Alu finished. Time required {} min\n'.format (round ((time.time() - start_time)/60),2))

######################### Section 2: Analysis of the data

# Define variables
path = "./"
output = os.listdir(path)
jacusa_out = [file for file in output if "Jacusa_" in file]

for file in jacusa_out:
    # Analysis
    ## Load data
    data = load_file(file)
    data["WTID"] = data['Region'].str.replace("*", ".") + '|' + data['Position'].astype(str) + '|' + data['wt'].str.replace('->', '')
    data["ADARID"] = data['Region'].str.replace("*", ".") + '|' + data['Position'].astype(str) + '|' + data['adarko'].str.replace('->', '')

    wt_dict = dict (zip(data.WTID, data.WTID))
    adar_dict = dict (zip(data.ADARID, data.ADARID))
    
    # Remove variants from WT present in ADARKO
    wt_clean = [res for res in wt_dict if res not in adar_dict and 'no change' not in res]
    adar_clean = [res for res in adar_dict if res not in wt_dict and 'no change' not  in res]
    
    # Remove variants present in dbSNP154
    wt_clean_snp = [ res for res in wt_clean if res not in dbSNP]
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
    
    support = file[-5]
    # Print output
    print ('\nStarting analysis on Jacusa2 output taking 6 samples (x3 WT_mock) & (x3 ADAR1KO_mock)\n')
    print ('\nStep 1: Load data (support ={}))\n'.format(support))
    print ('There are {} RES detected'.format(len(data)))
    print ('\nStep 2: Remove WT RES present in ADAR1KO\n')
    print ('There are {} RES in WT & {} in ADAR1KO'.format(len(wt_clean), len(adar_clean)))
    print ('\nStep 4: Remove RES present in dbSNP154\n')
    print ('There are {} RES in WT & {} in ADAR1KO'.format(len(wt_clean_snp), len(adar_clean_snp)))
    print ('{} WT RES and {} ADAR1KO RES were removed'.format(len (wt_clean) - len (wt_clean_snp), len(adar_clean)-len(adar_clean_snp)))
    print ('\nStep 5: Calculate % of RES present in REDIportal\n')
    print ('{} % of WT RES and {} % of ADAR1KO RES are in REDIportal'.format(round (len(wt_db)/len(A2G_T2C_wt)*100,5), round (len(adar_db)/len(A2G_T2C_adar)*100,5)))
    print ('\nStep 6: Calculate number of A > G & T > C RES\n')
    print ('There are {} A > G & {} T > C in WT sample'.format(len (A2G_wt), len(T2C_wt)))
    print ('\nStep 7: Calculate FDR G2A / (T2C+A2G)')
    print('There are {} G>A and the FDR is {}'.format(len(G2A), FDR))
    print ('\n\nFirst part of the analysis finished, moving to computing % in Alu regions\n\n')
    cont = 0
    for ID in A2G_T2C_wt:
        chrom = ID.split('|')[0]
        pos = ID.split('|')[1]
        if chrom in Alu:
            if pos in Alu[chrom]:
                cont += 1
    print ('There are {} Percentage of IDs in Alu repeats'.format(round (cont/len(A2G_T2C_wt)*100,2)))
  
print ("Analysis finished")


