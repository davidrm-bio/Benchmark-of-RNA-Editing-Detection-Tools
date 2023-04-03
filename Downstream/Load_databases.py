# Script to load databases

## Modules
import os
import ujson as json

## Functions

def load_database (path, file):
    '''Load JSON file into dictionary '''
    with open (os.path.join(path, file), 'r') as file:
        db = json.load (file)
    return db

def generate_data ():
    dbSNP = load_database ("~data/genome/", 'dbSNP109.json' )
    REDIportal = load_database ("~data/dbRNA-Editing", 'REDIportal.json')
    Alu = load_database ("~data/genome/", 'Alu_GRCh38.json')
    return dbSNP, REDIportal, Alu