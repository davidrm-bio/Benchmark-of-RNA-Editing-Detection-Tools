# Modules
import os
import pandas as pd
import numpy as np
import ujson as json
import warnings
import time
from load_databases import generate_data

warnings.filterwarnings("ignore")

# Functions
def load_file (filename):
    df = pd.read_csv(filename, delimiter= '2', skiprows=1) 
    return df

######################### Section 1: Load databases
start_time = time.time()

print ('\nLoading databases into memory, this will take time (~30 min)...')

dbSNP, REDIportal, Alu = generate_data() # Change files to dbSNP138.json, REDIportal_hg19.json and Alu_GRCh37.json

print('Loading of dbSNP138, REDIportal & Alu finished. Time required {} min\n'.format (round ((time.time() - start_time)/60),2))

######################### Section 2: Analysis of the data

try:
    print ('\n\n############################################')
    pval = float(input('Pval (Enter 0 to exit): '))
    print ('############################################\n\n')
except:
    print('ValueError setting to default pval = 0.5')
    pval = 0.5

## Define variables
path = "./"
output = os.listdir(path)
file = "RNA_editing.sites.txt"
cols = ['#Chromosome', 'Position', 'RefEd']
Adar_out = [file for file in output if "ADAR" in file]
WT_out = [file for file in output if "WT" in file]

while pval != 0:
    # Combine data and create ID dictionary
    ## WT
    wt_dict = dict()
    for sample in WT_out:
        data = load_file (os.path.join(sample, file))
        data_filt = data[data["P_edit"]>=pval]
        data_filt['RefEd'] = data_filt['Reference'] + data_filt['Alternative']
        data_filt['ID'] = data_filt[cols].apply(lambda row: '|'.join(row.values.astype(str)), axis=1)
        data_dict = dict(zip(data_filt.ID, data_filt.ID))
        for res in data_dict:
            if res in wt_dict:
                wt_dict[res] += 1
            else:
                wt_dict[res] = 1
    ## ADAR
    adar_dict = dict()
    for sample in Adar_out:
        data = load_file (os.path.join(sample, file))
        data_filt = data[data["P_edit"]>=pval]
        data_filt['RefEd'] = data_filt['Reference'] + data_filt['Alternative']
        data_filt['ID'] = data_filt[cols].apply(lambda row: '|'.join(row.values.astype(str)), axis=1)
        data_dict = dict(zip(data_filt.ID, data_filt.ID))
        for res in data_dict:
            if res in adar_dict:
                adar_dict[res] += 1
            else:
                adar_dict[res] = 1

    # Retain variants present in at least 2 clones
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
    
    # Correction of IDs
    wt_clean_corrected = list()
    for region in wt_clean:
        save = region.split('|')[0]
        if '.' in save:
            newID = region.replace('chr', '')
        else:
            newID = region
        wt_clean_corrected.append(newID)
    
    adar_clean_corrected = list()
    for region in adar_clean:
        save = region.split('|')[0]
        if '.' in save:
            newID = region.replace('chr', '')
        else:
            newID = region
        adar_clean_corrected.append(newID)

    # Remove possible SNPs
    wt_clean_snp = [res for res in wt_clean_corrected if res not in dbSNP]
    adar_clean_snp = [res for res in adar_clean_corrected if res not in dbSNP]

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

    # Print output
    print ('\nStarting analysis on RED-ML output taking 6 samples (x3 WT_mock) & (x3 ADAR1KO_mock)\n')
    print ('\nStep 1: Load data + Filtering (P-value >= {})\n'.format(pval))
    print ('There are {} RES in WT & {} RES in ADAR1KO'.format(len(wt_dict), len(adar_dict)))
    print ('\nStep 2: Retain RES present in at least 2 clones\n')
    print ('There are {} RES in WT & {} in ADAR1KO'.format(len(wt_clean_tmp), len(adar_clean_tmp)))
    print ('\nStep 3: Remove WT RES present in ADAR1KO\n')
    print ('There are {} RES in WT & {} in ADAR1KO'.format(len(wt_clean), len(adar_clean)))
    print ('\nStep 4: Remove RES present in dbSNP138\n')
    print ('There are {} RES in WT & {} in ADAR1KO'.format(len(wt_clean_snp), len(adar_clean_snp)))
    print ('{} WT RES and {} ADAR1KO RES were removed'.format(len (wt_clean) - len (wt_clean_snp), len(adar_clean)-len(adar_clean_snp)))
    print ('\nStep 5: Calculate % of RES present in REDIportal\n')
    print ('{} % of WT RES and {} % of ADAR1KO RES are in REDIportal'.format(round (len(wt_db)/len(A2G_T2C_wt)*100,5), round (len(adar_db)/len(A2G_T2C_adar)*100,5)))
    print ('\nStep 6: Calculate number of A > G & T > C RES\n')
    print ('There are {} A > G & {} T > C in WT sample'.format(len (A2G_wt), len(T2C_wt)))
    
    print ('\n\nFirst part of the analysis finished, moving to computing % in Alu regions\n\n')
    cont =0
    for ID in A2G_T2C_wt:
        chrom = ID.split('|')[0]
        pos = ID.split('|')[1]
        if chrom in Alu:
            if pos in Alu[chrom]:
                cont += 1

    print ('There are {} Percentage of IDs in Alu repeats'.format(round (cont/len(A2G_T2C_wt)*100,2)))
    
    try:
        print ('\n\n############################################')
        pval = float(input('Pval (Enter 0 to exit): '))
        print ('############################################\n\n')
    except:
        print('ValueError setting to default pval = 0.5')
        pval = 0.5

        
print ("Analysis finished")
