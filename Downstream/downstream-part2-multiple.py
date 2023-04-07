### Analysis considering clones part 2 ###
# Modules
import ujson as json
import pandas as pd
import numpy as np
from scipy.stats import sem
import os

# Functions 
def load_database (path, file):
    '''Load JSON file into dictionary '''
    with open (os.path.join(path, file), 'r') as file:
        db = json.load (file)
    return db

def generate_data ():
    #REDIportal = load_database ("/binf-isilon/rennie/gsn480/data/compare_tools/dataset1/dbRNA-Editing", 'REDIportal_hg19.json')
    #Alu = load_database ("/binf-isilon/rennie/gsn480/data/compare_tools/dataset1", 'Alu_GRCh37.json')
    
    REDIportal = load_database ("/binf-isilon/rennie/gsn480/data/compare_tools/dataset1/dbRNA-Editing", 'REDIportal.json')
    Alu = load_database ("/binf-isilon/rennie/gsn480/data/compare_tools/dataset1", 'Alu_GRCh38.json')
    return REDIportal, Alu

# Global variables
file = 'REDItools2_data.json' # Change accordingly ############################################################
aligners = ['bwa', 'hisat2', 'star'] # ['bwa', 'hisat2']
conditions = ['WT', 'ADAR1KO']
samples = ['clone1', 'clone2', 'clone3']
supports =  ['2','4','6','8','10'] # or ['0', '0.1'] or [2,4,6,8,10] or ['0.5','0.6','0.7','0.8','0.9']
measures = ['Variants', 'RES', 'FDR', 'REDIportal', 'Alu', 'List_of_RES', 'List_of_Vars']

# Load data
with open(file, 'r') as fp:
    data = json.load(fp)

print ('\n\nPre-Analysis to keep variants present in the 3 replicates')

# Pre-Analysis
data_clean = {align:{condition:{str (support):{'List_of_Vars': set(), 'List_of_RES': set()} for support in supports} 
                     for condition in conditions} for align in aligners}

# Clean the data keeping RES present in the 3 replicates - Pre-processing
for align in aligners:
    for condition in conditions:
        for support in supports:
            # Keep variants present in clone1 and clone2 and clone 3
            a = set (data[align][condition]['clone1'][support]['List_of_Vars'])
            b = set (data[align][condition]['clone2'][support]['List_of_Vars'])
            c = set (data[align][condition]['clone3'][support]['List_of_Vars'])
            common = a & b & c
            common_res = [res for res in common if 'AG' in res or 'TC' in res]
            data_clean[align][condition][support]['List_of_Vars'] = common
            data_clean[align][condition][support]['List_of_RES'] = common_res


print ('Finished... Proceed to the analysis')

# Analysis 
export_data = {align : {condition : {support: {measure:0 for measure in measures} 
                                     for support in supports} for condition in conditions} for align in aligners}

# Load data

REDIportal, Alu = generate_data()

with open ('REDItools2_multiple_Excel.csv', 'a') as excel: # Change accordingly ############################################################
    for align in aligners:
        for condition in conditions:
            for support in supports:
                #print ('Analysis ({} | {} | {})'.format(align, condition, support))
                data_res_tmp = list (data_clean[align][condition][support]['List_of_RES'])
                data_vars_tmp = list (data_clean[align][condition][support]['List_of_Vars'])

                # Check how many are in REDIportal
                data_db = [res for res in data_res_tmp if res in REDIportal]

                try:
                    db_percentage = round (len (data_db) / len(data_res_tmp) * 100, 3)
                except ZeroDivisionError:
                    db_percentage = 0

                # FDR
                data_GA = [var for var in data_vars_tmp if 'GA' in var]
                try:
                    data_fdr = len(data_GA) / len (data_res_tmp) * 100
                except:
                    data_fdr = 0

                # Count how many are in Alu regions
                cont =0
                for ID in data_res_tmp:
                    chrom = ID.split('|')[0]
                    pos = ID.split('|')[1]
                    if chrom in Alu:
                        if pos in Alu[chrom]:
                            cont += 1
                try:
                    alu_percentage = round (cont / len (data_res_tmp) * 100, 3)
                except ZeroDivisionError:
                    alu_percentage = 0


                # Save data
                export_data[align][condition][support]['Variants'] = len (data_vars_tmp)
                export_data[align][condition][support]['RES'] = len (data_res_tmp)
                export_data[align][condition][support]['REDIportal'] = len (data_db)
                export_data[align][condition][support]['Alu'] = cont
                export_data[align][condition][support]['FDR'] = data_fdr
                export_data[align][condition][support]['List_of_RES'] = data_res_tmp
                export_data[align][condition][support]['List_of_Vars'] = data_vars_tmp
                
                # Write for excel
                n_res = len (data_res_tmp)
                n_fdr = data_fdr
                n_rediportal = db_percentage
                n_alu = alu_percentage
                
                line = [condition, align.upper(), support, str(n_res), str(n_fdr), str(n_rediportal), str(n_alu)]
                
                # write into a file
                excel.write(','.join(line) + '\n')
    

print ('\n\nSaving data')

with open("REDItools2_multiple_data.json" , 'w') as fp: # Change accordingly  ############################################################
    json.dump(export_data, fp)
    
print ('Analysis finished')
            
            
            
