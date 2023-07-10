# Modules
import os
import ujson as json
import pandas as pd
from MainFunctions import  load_file, load_database
import re

# Global variables
path = "../JACUSA2"
aligners = ['bwa', 'hisat2', 'star']
conditions = ['WT', 'ADAR1KO']
samples = ['clone1', 'clone2', 'clone3']
categories = ['N_RES', 'N_REDIportal', 'N_Alu', 'SNPs', 'List_of_Res']
thresholds = ['2', '4', '6', '8', '10']
RES = ['AG', 'TC']
dbpath = '/binf-isilon/rennie/gsn480/data/compare_tools/dataset1/dbRNA-Editing'
hek = load_database(dbpath, 'HEK293T_hg38_clean.json')
rediportal38 = load_database(dbpath,'REDIportal.json')
# Dictionary for part 3
export_data = {align: {condition: {cutoff: { category:0 for category in categories} 
                                for cutoff in thresholds} 
            for condition in conditions}  
    for align in aligners}
# Analysis
for align in aligners:
    newpath = os.path.join(path, align)
    files = [file for file in os.listdir(newpath)]
    for file in files:
        # Load data
        data = load_file(os.path.join (newpath, file))
        # Retain only RES - Part 1 
        if align == 'hisat2':
            data["WTID"] = 'chr' + data['Region'].str.replace("*", ".") + '|' + data['Position'].astype(str) + '|' + data['wt'].str.replace('->', '')
            data["ADARID"] = 'chr' + data['Region'].str.replace("*", ".") + '|' + data['Position'].astype(str) + '|' + data['adarko'].str.replace('->', '')
        else: 
            data["WTID"] = data['Region'].str.replace("*", ".") + '|' + data['Position'].astype(str) + '|' + data['wt'].str.replace('->', '')
            data["ADARID"] = data['Region'].str.replace("*", ".") + '|' + data['Position'].astype(str) + '|' + data['adarko'].str.replace('->', '')
            
        wt  = [res for res in data.WTID.to_list() if res.split('|')[-1] in RES]
        adar  = [res for res in data.ADARID.to_list() if res.split('|')[-1] in RES]
        # Remove SNPs - Part 2
        def remove_SNPs (data):
            '''
            Providing a list, remove RES present in HEK dict
            '''
            clean = []
            cont = 0
            for res in data:
                if res in hek:
                    cont += 1
                else:
                    clean.append(res)
            return clean, cont
        wt_clean, wt_snps = remove_SNPs(wt)
        adar_clean, adar_snps = remove_SNPs(adar)
        # Summarise data - part 3
        if align == 'hisat2':
            cutoff = re.findall(r'\d+', file)[1]
        else:
            cutoff = re.findall(r'\d+', file)[0]
        # RES in REDIportal
        wt_db = [res for res in wt_clean if res in rediportal38]
        adar_db = [res for res in adar_clean if res in rediportal38]
        # Save WT
        export_data[align]['WT'][cutoff]['N_RES'] = len (wt_clean)
        export_data[align]['WT'][cutoff]['N_REDIportal'] = len (wt_db)
        export_data[align]['WT'][cutoff]['SNPs'] = wt_snps
        export_data[align]['WT'][cutoff]['List_of_Res'] = wt_clean
        # Save ADAR1KO
        export_data[align]['ADAR1KO'][cutoff]['N_RES'] = len (adar_clean)
        export_data[align]['ADAR1KO'][cutoff]['N_REDIportal'] = len (adar_db)
        export_data[align]['ADAR1KO'][cutoff]['SNPs'] = adar_snps
        export_data[align]['ADAR1KO'][cutoff]['List_of_Res'] = adar_clean
        
# Save  data in JSON
file_to_save = 'Data_JACUSA2.json'
with open (os.path.join ("../", file_to_save), 'w') as out:
    json.dump(export_data, out)


        
            
                    
        
        
