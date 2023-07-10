# Modules
import os
import ujson as json
from MainFunctions import load_database, stats


# Global Variables
path = "../BCFTools"
aligners = ['bwa', 'hisat2', 'star']
conditions = ['WT', 'ADAR1KO']
samples = ['clone1', 'clone2', 'clone3']
RES = ['AG', 'TC']
output_dir = 'ResultsFiles'
categories = ['N_RES', 'N_REDIportal', 'N_Alu', 'SNPs', 'List_of_Res']
thresholds = ['0', '0.01']

dbpath = '/binf-isilon/rennie/gsn480/data/compare_tools/dataset1/dbRNA-Editing'
hek = load_database(dbpath,'HEK293T_hg38_clean.json')
rediportal38 = load_database(dbpath,'REDIportal.json')

'''
# Part 1 - Convert VCF to Bed
for align in aligners:
    newpath = os.path.join(path, align)
    if  os.path.exists(output_dir) == False:
        os.system('mkdir {}'.format(os.path.join(newpath, output_dir)))
    files = [file for file in os.listdir(newpath) if 'vcf' in file]
    for file in files:
        newfile = file.split('filt')[0] + file.split('sortrmdup_2.')[1] 
        newfile = newfile.replace('.vcf', '.bed')
        os.system('vcf2bed < {} > {}'.format(os.path.join (newpath,file), os.path.join(newpath, output_dir,newfile))) # Convert tool

# Part  2 - Keep RES ; 
for align in aligners:
    newpath = os.path.join(path, align, output_dir)
    files = [file for file in os.listdir(newpath) if 'bed' in file]
    for file in files: 
        newfile = file.replace('.bed', '_processed.bed')
        with open (os.path.join(newpath, newfile), 'w') as out:
            with open (os.path.join (newpath, file), 'r') as Read:
                for lines in Read:
                    line = lines.strip().split('\t')
                    res = line[5] + line[6]
                    if res in RES:
                        if align == 'hisat2':
                            line[0] = 'chr' + line[0]
                        out.write ('\t'.join(line) + '\n')

# Part 3 - Remove HEK293T cells SNPs
for align in aligners:
    newpath = os.path.join(path, align, output_dir)
    files = [file for file in os.listdir(newpath) if '_processed.bed' in file]
    for file in files:
        newfile = file.replace('_processed.bed', '_processed_HEK.bed')
        cont = 0
        with open (os.path.join(newpath, newfile), 'w') as out:
            with open (os.path.join (newpath, file), 'r') as Read:
                for lines in Read:
                    line = lines.strip().split('\t')
                    ids = line[0] + '|' + line[2] + '|' + line [5] + line [6]
                    if ids in hek:
                        cont +=1
                    else:
                        out.write (lines)
                out.write ('# A total of {} RES were removed (Present in HEK293T)'.format(cont))
'''

# Part 4 - Summarise data
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
                if cutoff == '0':
                    cutoff_file = '_0_processed'
                else:
                    cutoff_file = '0_0_1_processed'
                file = [files for files in os.listdir(newpath) if condition in files and clone in files and 'HEK' in files and cutoff_file in files][0]
                res, db, snp, res_list = stats(os.path.join  (newpath, file), support= float(cutoff), tool = 'BCFTools', rediportal=rediportal38)
                # Save data
                data[align][condition][clone][cutoff]['N_RES'] = res
                data[align][condition][clone][cutoff]['N_REDIportal'] = db
                data[align][condition][clone][cutoff]['SNPs'] = snp
                data[align][condition][clone][cutoff]['List_of_Res'] = res_list
# Save  data in JSON
file_to_save = 'Data_BCFTools.json'
with open (os.path.join ("../", file_to_save), 'w') as out:
    json.dump(data, out)
