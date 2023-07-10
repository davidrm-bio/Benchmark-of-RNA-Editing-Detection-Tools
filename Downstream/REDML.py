# Modules
import os
import ujson as json
from MainFunctions import load_database, stats

# Global Variables
path = "../REDML"
aligners = ['bwa', 'hisat2', 'star']
conditions = ['WT', 'ADAR1KO']
samples = ['clone1', 'clone2', 'clone3']
Pval, RES = 0.5, ['AG', 'TC']
output_file = "RNA_editing.sites.txt"
output_dir = 'ResultsFiles'
categories = ['N_RES', 'N_REDIportal', 'N_Alu', 'SNPs', 'List_of_Res']
thresholds = ['0.5', '0.6', '0.7', '0.8', '0.9']

dbpath = '/binf-isilon/rennie/gsn480/data/compare_tools/dataset1/dbRNA-Editing'
hek = load_database(dbpath ,'HEK293T_hg19_clean.json')
rediportal37 = load_database(dbpath,'REDIportal_hg19.json')

'''
# Part 1 - Keep RES ; Supporting reads >= 2 
for align in aligners:
    newpath  = os.path.join(path, align)
    if  os.path.exists(output_dir) == False:
        os.system('mkdir {}'.format(os.path.join(newpath, output_dir)))
    for file in os.listdir(newpath):
        if file == output_dir:
            continue
        newfile = file.replace('_output', '_processed.output')
        with open (os.path.join(newpath,output_dir, newfile), 'w') as out:
            with open (os.path.join (newpath, file, output_file), 'r') as Read:
                cont = 0
                for lines in Read:
                    cont += 1 
                    if cont == 1:
                        continue
                    line = lines.strip().split('\t')
                    res = line[3] + line[5]
                    if res in RES and float (line[7]) >= Pval:
                        out.write ('\t'.join(line) + '\n')

# Part 2 - Remove HEK293T cells SNPs
for align in aligners:
    newpath = os.path.join(path, align, output_dir)
    files = [file for file in os.listdir(newpath) if 'processed' in file]
    for file in files:
        newfile = file.replace('_processed.output', '_processed_HEK.output')
        cont = 0
        with open (os.path.join(newpath, newfile), 'w') as out:
            with open (os.path.join (newpath, file), 'r') as Read:
                for lines in Read:
                    line = lines.strip().split('\t')
                    ids = line[0] + '|' + line[1] + '|' + line [3] + line [5]
                    if ids in hek:
                        cont +=1
                    else:
                        out.write (lines)
                out.write ('# A total of {} RES were removed (Present in HEK293T)'.format(cont))
'''

# Part 3 - Summarise data
# Dictionary to save the data
data = {align: {condition: {clone: {cutoff: { category:0 for category in categories} 
                                for cutoff in thresholds} 
                        for clone in samples} 
            for condition in conditions}  
    for align in aligners}
# Calculate data
for align in aligners:
    newpath = os.path.join(path, align, output_dir)
    for condition in conditions:
        for clone in samples:
            for cutoff in thresholds:
                file = [files for files in os.listdir(newpath) if condition in files and clone in files and 'HEK' in files][0]
                res, db, snp, res_list = stats(os.path.join  (newpath, file), support= float(cutoff), tool = 'REDML', rediportal=rediportal37)
                # Save data
                data[align][condition][clone][cutoff]['N_RES'] = res
                data[align][condition][clone][cutoff]['N_REDIportal'] = db
                data[align][condition][clone][cutoff]['SNPs'] = snp
                data[align][condition][clone][cutoff]['List_of_Res'] = res_list
# Save  data in JSON
file_to_save = 'Data_REDML.json'
with open (os.path.join ("../", file_to_save), 'w') as out:
    json.dump(data, out)
