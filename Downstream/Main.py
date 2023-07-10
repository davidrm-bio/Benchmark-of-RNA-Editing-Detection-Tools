# Modules
import  os
import ujson as json
from MainFunctions import *

# Global Variables
path = "../"
aligners = ['bwa', 'hisat2', 'star']
conditions = ['WT', 'ADAR1KO']
samples = ['clone1', 'clone2', 'clone3']
tools = ['BCFTools', 'SPRINT', 'REDItools2', 'REDML']
categories = ['N_RES', 'N_REDIportal', 'N_Alu', 'SNPs', 'List_of_Res']
dbpath = '/binf-isilon/rennie/gsn480/data/compare_tools/dataset1/dbRNA-Editing'

# Run Processing Scripts - Single analysis
os.system('python REDItools2.py')  
os.system('python SPRINT.py')          
os.system('python REDML.py')
os.system('python BCFtools.py') 
os.system('python JACUSA2.py')

# Update with Alu counts
os.system('python  Alu.py')

# Individual analysis
os.system('python Individual-Analysis.py')

# Multiple Analaysis
os.system('python Re-Analysis-Multiple.py')
os.system('python Multiple-Analysis.py')


