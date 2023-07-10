# Modules
import pandas as pd 
import os
import ujson as json

# Functions

def hek_SNPs ():
    '''
    Return a dictionary with known A > G or T > C SNPs for HEK293T
    Format of input BED file
    chr  pos0    pos1     ID (Chr|pos1|RES)
    '''
    with open ('/binf-isilon/rennie/gsn480/data/compare_tools/dataset1/dbRNA-Editing/HEK293T_hg38_clean.bed', 'r') as Read:
        hek = [lines.strip().split('\t')[-1] for lines in Read]      

    hek = dict (zip(hek, hek))
    return hek


# Load and create ID on BED file 
with open ('HEK293T_hg38_clean.bed',  'w') as out:
    with open ('HEK293T_hg38.bed', 'r') as Read:
        for lines in Read:
            line = lines.split('\t')
            chrom = line[0]
            pos = line[2]
            ref = line[3][-2]
            alt = line[3][-1]
            ID = chrom + '|' + pos + '|' + ref + alt
            newline = [chrom, str (int (pos) - 1), pos, ID]
            out.write('\t'.join(newline) + '\n')
            
# Generate dict for HEK293T
Hek = hek_SNPs()

# Save as Json file for future analysis
with open("'HEK293T_hg38_clean.json'", "w") as fp:
    json.dump(Hek , fp) 
