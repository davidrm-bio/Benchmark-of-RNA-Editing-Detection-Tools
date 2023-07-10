# Modules
import ujson as json
from MainFunctions import  load_database, count_Alu, save2JSON


# Global variables
path = '../'
dbpath = '/binf-isilon/rennie/gsn480/data/compare_tools/dataset1/dbRNA-Editing'
print ('Loading Alu GCh38...')
alu38 = load_database(dbpath,'Alu_GRCh38.json')
print ('Finished')
aligners = [['bwa', 'hisat2', 'star'],  # For the rest
            ['bwa', 'hisat2']] # For SPRINT
samples = ['clone1', 'clone2', 'clone3']
conditions = ['WT', 'ADAR1KO']
thresholds = [['2', '4', '6', '8', '10'],  # For the rest
             ['0.5','0.6', '0.7', '0.8', '0.9'], # For REDML
             ['0', '0.01']] # For BCFtools - Correction needed 0.01 -> 0.1

# REDItools2
reditools = load_database(path, 'Data_REDItool2.json')
for align in aligners[0]:
    for condition in conditions:
        for clone in samples:
            for cutoff in thresholds[0]:
                cont_reditools = count_Alu (reditools[align][condition][clone][cutoff]['List_of_Res'], alu38)
                reditools[align][condition][clone][cutoff]['N_Alu'] = cont_reditools

save2JSON (path, 'Data_REDItool2.json', reditools)

print ('REDITOOLS2 DONE')

# JACUSA2
jacusa = load_database(path, 'Data_JACUSA2.json')
for align in aligners[0]:
    for condition in conditions:
        for cutoff in thresholds[0]:
            cont_jacusa = count_Alu (jacusa[align][condition][cutoff]['List_of_Res'], alu38)
            jacusa[align][condition][cutoff]['N_Alu'] = cont_jacusa
            
save2JSON (path, 'Data_JACUSA2.json', jacusa)

print ('JACUSA2 DONE')

# SPRINT
sprint = load_database(path, 'Data_SPRINT.json')
for align in aligners[1]:
    for condition in conditions:
        for clone in samples:
            for cutoff in thresholds[0]:
                cont_sprint = count_Alu (sprint[align][condition][clone][cutoff]['List_of_Res'], alu38)
                sprint[align][condition][clone][cutoff]['N_Alu'] = cont_sprint

save2JSON (path, 'Data_SPRINT.json', sprint)

print ('SPRINT DONE')
# BCFTools
bcftools = load_database(path, 'Data_BCFTools.json')
for align in aligners[0]:
    for condition in conditions:
        for clone in samples:
            for cutoff in thresholds[2]:
                cont_bcftools = count_Alu (bcftools[align][condition][clone][cutoff]['List_of_Res'], alu38)
                bcftools[align][condition][clone][cutoff]['N_Alu'] = cont_bcftools

save2JSON (path, 'Data_BCFTools.json', bcftools)

print ('BCFTools DONE')


print ('Loading Alu GCh37...')
alu37 = load_database(dbpath,'Alu_GRCh37.json')
print ('Finished...')

# REDML
redml = load_database(path, 'Data_REDML.json')
for align in aligners[0]:
    for condition in conditions:
        for clone in samples:
            for cutoff in thresholds[1]:
                cont_redml = count_Alu (redml[align][condition][clone][cutoff]['List_of_Res'], alu37)
                redml[align][condition][clone][cutoff]['N_Alu'] = cont_redml

save2JSON (path, 'Data_REDML.json', redml)



