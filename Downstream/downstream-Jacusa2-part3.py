### Individual analysis part 2 ###
# Modules
import ujson as json
import pandas as pd
import numpy as np
from scipy.stats import sem


# Global variables
file = 'JACUSA2_data.json' # Change accordingly
aligners = ['bwa', 'hisat2', 'star'] # ['bwa', 'hisat2']
conditions = ['wt', 'adarko']
supports =  [2,4,6,8,10]  
measures = ['Variants', 'RES', 'FDR', 'REDIportal', 'Alu']


# Load data
with open(file, 'r') as fp:
    data = json.load(fp)
    
# Stats

with open ('JACUSA2_Excel.csv', 'a') as excel: # Change accordingly
    for align in aligners:
        for condition in conditions:
            for support in supports:
                support = str(support)
                variants = data[align][condition][support]['Variants']
                RES = data[align][condition][support]['RES']
                FDR = data[align][condition][support]['FDR']
                REDIportal = data[align][condition][support]['REDIportal']
                Alu = data[align][condition][support]['Alu']

                REDIportal_perc = REDIportal / RES * 100
                Alu_perc =  Alu / RES * 100
                
                # Write for excel
                n_res = np.round (RES, 0)
                n_fdr = np.round (FDR, 5)
                n_rediportal = np.round (REDIportal_perc, 2)
                n_alu = np.round (Alu_perc,2 )
                
                line = [condition, align.upper(), support, str(n_res), str(n_fdr), str(n_rediportal), str(n_alu)]
                
                # write into a file
                excel.write(','.join(line) + '\n')
    
    
    
