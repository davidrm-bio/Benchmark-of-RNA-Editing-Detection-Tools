### Analysis for JACUSA2 ###
# RES reported in WT and ADAR1KO samples are present in the 3 replicates (robust function from JACUSA2helper was used)
# Modules
import os
import pandas as pd
import numpy as np
import ujson as json
import warnings
import time
from MyFunctions import load_file
from load_databases import generate_data

# Load databases
start_time = time.time()
print ('\nLoading databases into memory, this will take time (~20 min)...')
dbSNP, REDIportal, Alu = generate_data()
print('Loading of dbSNP154, REDIportal & Alu finished. Time required {} min\n'.format (round ((time.time() - start_time)/60),2))

# Global variables
path = "./"
aligners = ['bwa', 'hisat2', 'star']
conditions = ['wt', 'adarko']
supports = [2, 4, 6, 8, 10]
measures = ['Variants', 'RES', 'FDR', 'REDIportal', 'Alu', 'List_of_RES', 'List_of_Vars']

export_data = {align : {condition : {support: {measure:0 for measure in measures} for support in supports} 
                        for condition in conditions} for align in aligners}

# Analysis 
for align in aligners:
    files = os.listdir(path) # ./ + bwa or hisat2 or star + /
    
    for support in supports:
        # Load File
        file = [f for f in files if  align in f and str(support) +'.tab' in f][0] 

        # Load data & filter by supporting reads
        data = load_file(file)
        data["WTID"] = data['Region'].str.replace("*", ".") + '|' + data['Position'].astype(str) + '|' + data['wt'].str.replace('->', '')
        data["ADARID"] = data['Region'].str.replace("*", ".") + '|' + data['Position'].astype(str) + '|' + data['adarko'].str.replace('->', '')
        
        # Convert to list
        data_wt = data.WTID.to_list()
        data_adar = data.ADARID.to_list()
        
        # Correction for Hisat2
        if align == 'hisat2': 
            data_wt_corrected_tmp = []
            for var in data_wt:
                if '.' in var: 
                    newvar = var
                else:
                    newvar = 'chr' + var
                data_wt_corrected_tmp.append(newvar)
            
            data_adar_corrected_tmp = []
            for var in data_adar:
                if '.' in var: 
                    newvar = var
                else:
                    newvar = 'chr' + var
                data_adar_corrected_tmp.append(newvar)
            
            data_wt = data_wt_corrected_tmp
            data_adar = data_adar_corrected_tmp
            

        # Remove Variants in dbSNP
        data_wt_ID_clean = [res for res in data_wt if res not in dbSNP]
        data_adar_ID_clean = [res for res in data_adar if res not in dbSNP]

        # Get number of RES
        data_wt_res = [res for res in data_wt_ID_clean if 'AG' in res or 'TC' in res]
        data_adar_res = [res for res in data_adar_ID_clean if 'AG' in res or 'TC' in res]

        # Check how many are in REDIportal
        data_wt_db = [res for res in data_wt_res if res in REDIportal]
        data_adar_db = [res for res in data_adar_res if res in REDIportal]
        
        try:
            db_wt_percentage = round (len (data_wt_db) / len(data_wt_res) * 100, 3)
            db_adar_percentage = round (len (data_adar_db) / len(data_adar_res) * 100, 3)
        except ZeroDivisionError:
            db_wt_percentage = 0
            db_adar_percentage = 0

        # FDR
        data_wt_GA = [var for var in data_wt_ID_clean if 'GA' in var]
        data_adar_GA = [var for var in data_adar_ID_clean if 'GA' in var]
        
        try:
            data_wt_fdr = len(data_wt_GA) / len (data_wt_res) * 100
            data_adar_fdr = len(data_adar_GA) / len (data_adar_res) * 100
        except:
            data_wt_fdr = 0
            data_adar_fdr = 0

        # Count how many are in Alu regions
        cont_wt =0
        for ID in data_wt_res:
            chrom = ID.split('|')[0]
            pos = ID.split('|')[1]
            if chrom in Alu:
                if pos in Alu[chrom]:
                    cont_wt += 1
        try:
            alu_wt_percentage = round (cont_wt / len (data_wt_res) * 100, 3)
        except ZeroDivisionError:
            alu_wt_percentage = 0
        
        cont_adar =0
        for ID in data_adar_res:
            chrom = ID.split('|')[0]
            pos = ID.split('|')[1]
            if chrom in Alu:
                if pos in Alu[chrom]:
                    cont_adar += 1
        try:
            alu_adar_percentage = round (cont_adar / len (data_adar_res) * 100, 3)
        except ZeroDivisionError:
            alu_adar_percentage = 0        
        
        # Save data
        export_data[align]['wt'][support]['Variants'] = len (data_wt_ID_clean)
        export_data[align]['wt'][support]['RES'] = len (data_wt_res)
        export_data[align]['wt'][support]['REDIportal'] = len (data_wt_db)
        export_data[align]['wt'][support]['Alu'] = cont_wt
        export_data[align]['wt'][support]['FDR'] = data_wt_fdr
        export_data[align]['wt'][support]['List_of_RES'] = data_wt_res
        export_data[align]['wt'][support]['List_of_Vars'] = data_wt_ID_clean

        export_data[align]['adarko'][support]['Variants'] = len (data_adar_ID_clean)
        export_data[align]['adarko'][support]['RES'] = len (data_adar_res)
        export_data[align]['adarko'][support]['REDIportal'] = len (data_adar_db)
        export_data[align]['adarko'][support]['Alu'] = cont_adar
        export_data[align]['adarko'][support]['FDR'] = data_adar_fdr
        export_data[align]['adarko'][support]['List_of_RES'] = data_adar_res
        export_data[align]['adarko'][support]['List_of_Vars'] = data_adar_ID_clean
                
print ('\n\nSaving data')

with open("JACUSA2_data.json" , 'w') as fp:
    json.dump(export_data, fp)
    
print ('Analysis finished')
