# Modules
import pandas as pd 
import os
import ujson as json

# Functions

def load_file (filename):
    '''Load tab-delimited files into a DataFrame'''
    df = pd.read_csv(filename, delimiter="\t") 
    return df

def preprocessing(db):
    '''Generate ID list of REDIportal'''
    db['RefEd'] = db['Ref'] + db['Ed']
    cols = ['Region', 'Position', 'RefEd']
    db['ID'] = db[cols].apply(lambda row: '|'.join(row.values.astype(str)), axis=1)
    db_id = db['ID'].to_list()
    db_id_dict = dict (zip(db_id, db_id))
    return db_id_dict


# Load data
REDIportal = load_file("REDIportal_db_GRCh38.txt")

# Process data to generate a dictionary with the format {ID: SNP ID}
# where ID is Chr|Pos|RefAlt

REDIportal_clean = preprocessing (REDIportal)

# Save as Json file for future analysis
with open("REDIportal.json", "w") as fp:
    json.dump(REDIportal_clean , fp) 
