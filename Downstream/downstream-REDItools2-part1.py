### Individual analysis for REDItools2 ###
# Modules
import os
import pandas as pd
import numpy as np
import ujson as json
import warnings
import time
from MyFunctions import load_file
from load_databases import generate_data

warnings.filterwarnings("ignore")

# Load databases
start_time = time.time()
print ('\nLoading databases into memory, this will take time (~20 min)...')
dbSNP, REDIportal, Alu = generate_data()
print('Loading of dbSNP154, REDIportal & Alu finished. Time required {} min\n'.format (round ((time.time() - start_time)/60),2))


# Global variables
path = "./"
aligners = ['bwa', 'hisat2', 'star']
conditions = ['WT', 'ADAR1KO']
samples = ['clone1', 'clone2', 'clone3']
supports = [2,4,6,8,10]
measures = ['Variants', 'RES', 'FDR', 'REDIportal', 'Alu', 'List_of_RES', 'List_of_Vars']
cols = ['Region', 'Position', 'AllSubs']

export_data = {align : {condition : {sample : {support: 
                                               {measure:0 for measure in measures} for support in supports} 
                                     for sample in samples} for condition in conditions} for align in aligners}


# Analysis 
for align in aligners:
    files = os.listdir(path + align) # ./ + bwa or hisat2 or star + /
    newpath = os.path.join (path + align)
    for condition in conditions:
        for sample in samples:
            for support in supports:
                # Load File
                file = [f for f in files if condition in f and sample in f and 'output' in f][0] 
                
                # Load data & filter by supporting reads
                data = load_file (os.path.join(newpath, file))
                data_filt = data[data["Coverage-q30"]>=support]
                data_filt = data_filt[data_filt["Frequency"]>=0.1]
                data_filt['ID'] = data_filt[cols[0]].astype(str) + '|' + data_filt[cols[1]].astype(str) + '|' + data_filt[cols[2]]
                
                # Convert to lis
                data_ID = data_filt.ID.to_list()
                
                # Correction; Some positions have multiple Variants
                data_ID_corrected = []
                for ID in data_ID:
                    change = ID.split('|')[2]
                    if len(change) > 2: 
                        new_ID = ID.split('|')[0] + '|' + ID.split('|')[1] + '|' + ID.split('|')[2][0:2]
                    else:
                        new_ID = ID
                    data_ID_corrected.append(new_ID)
                        
                # Correction for Hisat2
                if align == 'hisat2': 
                    data_ID_corrected_tmp = []
                    for var in data_ID:
                        if '.' in var: 
                            newvar = var
                        else:
                            newvar = 'chr' + var
                        data_ID_corrected_tmp.append(newvar)
                    data_ID_corrected = data_ID_corrected_tmp
                
                # Remove Variants in dbSNP
                data_ID_clean = [res for res in data_ID_corrected if res not in dbSNP]
                
                # Get number of RES
                data_res = [res for res in data_ID_clean if 'AG' in res or 'TC' in res]
                
                # Check how many are in REDIportal
                data_db = [res for res in data_res if res in REDIportal]
                try:
                    db_percentage = round (len (data_db) / len(data_res) * 100, 3)
                except ZeroDivisionError:
                    db_percentage = 0
                
                # FDR
                data_GA = [var for var in data_ID_clean if 'GA' in var]
                try:
                    data_fdr = len(data_GA) / len (data_res) * 100
                except:
                    data_fdr = 0
                
                # Count how many are in Alu regions
                cont =0
                for ID in data_res:
                    chrom = ID.split('|')[0]
                    pos = ID.split('|')[1]
                    if chrom in Alu:
                        if pos in Alu[chrom]:
                            cont += 1
                try:
                    alu_percentage = round (cont / len (data_res) * 100, 3)
                except ZeroDivisionError:
                    alu_percentage = 0
                
                '''
                # Print out results 
                print ('\n\n##################################################')
                print ('Analysis in REDItools2 ({} | {} | {} | {})'.format (align, condition, sample, support))
                print ('Number of Variants: {}'.format(len(data_ID_clean)))
                print ('Number of RES: {}'.format(len(data_res)))
                print ('Number of RES in REDIportal: {} ({} %)'.format(len(data_db), db_percentage))
                print ('Number of G > A: {}'.format(len(data_GA)))
                print ('FDR: {}'.format(data_fdr))
                print ('Number of RES in Alu: {} ({} %)'.format(cont, alu_percentage))
                print ('\n\n##################################################')
                '''
                # Save data
                export_data[align][condition][sample][support]['Variants'] = len (data_ID_clean)
                export_data[align][condition][sample][support]['RES'] = len (data_res)
                export_data[align][condition][sample][support]['REDIportal'] = len (data_db)
                export_data[align][condition][sample][support]['Alu'] = cont
                export_data[align][condition][sample][support]['FDR'] = data_fdr
                export_data[align][condition][sample][support]['List_of_RES'] = data_res
                export_data[align][condition][sample][support]['List_of_Vars'] = data_ID_clean
                
print ('\n\nSaving data')

with open("REDItools2_data.json" , 'w') as fp:
    json.dump(export_data, fp)
    
print ('Analysis finished')
