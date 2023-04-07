### Individual analysis part 2 ###
# Modules
import ujson as json
import pandas as pd
import numpy as np
from scipy.stats import sem


# Global variables
file = 'REDItools2_data.json' # Change accordingly
aligners = ['bwa', 'hisat2', 'star'] # ['bwa', 'hisat2']
conditions = ['WT', 'ADAR1KO']
samples = ['clone1', 'clone2', 'clone3']
supports = [2,4,6,8,10]  # [2,4,6,8,10] # or [0.5,0.6,0.7,0.8,0.9]   or ['0', '0.1']
measures = ['Variants', 'RES', 'FDR', 'REDIportal', 'Alu']


# Load data
with open(file, 'r') as fp:
    data = json.load(fp)
    
# Stats

with open ('REDItools2_Excel.csv', 'a') as excel: # Change accordingly
    for align in aligners:
        for condition in conditions:
            for support in supports:
                support = str(support)
                variants, RES, FDR, REDIportal, Alu = [], [], [], [], []
                for sample in samples:
                    variants.append(data[align][condition][sample][support]['Variants'])
                    RES.append(data[align][condition][sample][support]['RES'])
                    FDR.append(data[align][condition][sample][support]['FDR'])
                    REDIportal.append(data[align][condition][sample][support]['REDIportal'])
                    Alu.append(data[align][condition][sample][support]['Alu'])

                if len(variants) > 3 or len(RES) > 3 or len (FDR) > 3 or len (REDIportal) > 3 or len (Alu) > 3:
                    raise ValueError('More than 3 values detected')

                REDIportal_perc = np.array (REDIportal) / np.array (RES) * 100
                Alu_perc =  np.array (Alu) / np.array(RES) * 100

                # Print results
                '''
                print ('\n\n##################################################')
                print ('Analysis in REDItools2 ({} | {} | {})'.format (align, condition, support))
                print ('Average number of Variants: {} (Std {})'.format(np.round (np.mean(variants), 0), np.round (sem(variants), 3)))
                print ('Average number of RES: {} (Std {})'.format(np.round (np.mean(RES), 0), np.round (sem(RES), 3)))
                print ('Average number of RES in REDIportal: {} % (Std {})'.format(np.round (np.mean(REDIportal_perc), 2), np.round (sem(REDIportal_perc), 5)))
                print ('FDR: {} (Std {})'.format(np.round (np.mean(FDR), 5), np.round (sem(FDR), 5)))
                print ('Number of RES in Alu: {} % (Std {})'.format(np.round (np.mean(Alu_perc), 2), np.round (sem(Alu_perc), 5)))
                print ('\n\n##################################################')
                '''
                mean_res = np.round(np.mean(RES), 0)
                sem_res = np.round (sem(RES), 3)
                mean_fdr = np.round(np.mean(FDR), 5)
                mean_rediportal = np.round(np.mean(REDIportal_perc), 2)
                mean_alu = np.round(np.mean(Alu_perc), 2)
                
                line = [condition, align.upper(), support, str(mean_res), str(sem_res), str(mean_fdr), str(mean_rediportal), str(mean_alu)]
                
                # write into a file
                excel.write(','.join(line) + '\n')
    
    
