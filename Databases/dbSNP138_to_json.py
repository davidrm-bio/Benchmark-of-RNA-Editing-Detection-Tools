# Modules
import os
import pandas as pd
import ujson as json

# Functions
def load_file (filename):
    '''Load tab-delimited files into a DataFrame'''
    df = pd.read_csv(filename, delimiter="\t",header=None) 
    return df

# Load data
dbSNP = load_file("snp138.txt")

# Process data to generate a dictionary with the format {ID: SNP ID}
# where ID is Chr|Pos|RefAlt

dbSNP_filt = dbSNP[dbSNP[10] == 'genomic']
dbSNP_filt = dbSNP_filt[dbSNP_filt[11] == 'single']
dbSNP_filt["ID"] = dbSNP_filt[1].astype(str) + '|' +  dbSNP[3].astype(str) + '|' + dbSNP[9].str.replace("/", "")

dbSNP_temp = dbSNP_filt[['ID', 4]]
dbSNP_temp = dbSNP_temp.drop_duplicates('ID')

dbSNP_dict = dbSNP_temp.set_index('ID').to_dict('index')
dbSNP_clean = dict()
for key in dbSNP_dict: 
    dbSNP_clean[key] = dbSNP_dict[key][4]

# Save as Json file for future analysis
with open("dbSNP138.json", "w") as fp:
    json.dump(dbSNP_clean , fp) 
