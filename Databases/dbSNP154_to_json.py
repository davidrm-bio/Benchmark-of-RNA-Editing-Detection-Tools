# Modules
import os
import pandas as pd
import ujson as json

# Functions
def load_file (filename):
    '''Load tab-delimited files into a DataFrame'''
    df = pd.read_csv(filename, delimiter="\t", names=["chr", "pos", "ID_SNP", "Ref", "Alt", "V6", "V7", "Info"]) 
    return df

# Load data
dbSNP = load_file("dbSNP109_concat_snps.vcf")

# Process data to generate a dictionary with the format {ID: SNP ID}
# where ID is Chr|Pos|RefAlt

dbSNP["ID"] = 'chr' + dbSNP['chr'].astype(str) + '|' +  dbSNP['pos'].astype(str) + '|' + dbSNP['Ref']
dbSNP_temp = dbSNP[['ID', 'Alt', 'ID_SNP']]
dbSNP_dict = dbSNP_temp.set_index('ID').to_dict('index')

dbSNP_clean = dict()
for key in dbSNP_dict:    
    if len (dbSNP_dict[key]['Alt']) > 1:
        for base in dbSNP_dict[key]['Alt'].split(','):
            dbSNP_clean[key + base] = dbSNP_dict[key]['ID_SNP']
    else:
        dbSNP_clean[key + dbSNP_dict[key]['Alt']] = dbSNP_dict[key]['ID_SNP']

# Save as Json file for future analysis
with open("dbSNP109.json", "w") as fp:
    json.dump(dbSNP_clean , fp) 