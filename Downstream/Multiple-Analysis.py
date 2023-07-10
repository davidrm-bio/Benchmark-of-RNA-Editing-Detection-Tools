# Modules 
import ujson as json
from MainFunctions import  *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats as sem
from more_itertools import locate

# Global variables
path = '..'
aligners = [['bwa', 'hisat2', 'star'],  # For the rest
            ['bwa', 'hisat2']] # For SPRINT
conditions = ['WT', 'ADAR1KO']

thresholds = [['2', '4', '6', '8', '10'],  # For the rest
             ['0.5','0.6', '0.7', '0.8', '0.9'], # For REDML
             ['0', '0.1']] # For BCFtools
categories = ['N_res', 'SEM_res', 'N_db', 'SEM_db', 'N_Alu', 'SEM_Alu']
tools = ['BCFTools', 'RED-ML', 'SPRINT', 'REDItools2']

# Load Data
reditools = load_database(path, 'Data_REDItool2-Multiple.json')
redml = load_database(path, 'Data_REDML-Multiple.json')
sprint = load_database(path, 'Data_SPRINT-Multiple.json')
bcftools = load_database(path, 'Data_BCFTools-Multiple.json')
jacusa = load_database(path, 'Data_JACUSA2.json')

# Generate Tables
save_to_cs_multiple('REDItools2-Multiple_Table.csv', reditools, aligners[0], conditions, thresholds[0])
save_to_cs_multiple('SPRINT-Multiple_Table.csv', sprint, aligners[1], conditions, thresholds[0])
save_to_cs_multiple('REDML-Multiple_Table.csv', redml, aligners[0], conditions, thresholds[1], thresholdname='Detection Threshold')
save_to_cs_multiple('BCFTools-Multiple_Table.csv', bcftools, aligners[0], conditions, thresholds[2], thresholdname = 'Minor Allele Frequency')
save_to_csv_multiple('JACUSA2_Table.csv', jacusa, aligners[0], conditions, thresholds[0])


# Data for plotting
# REDItools2
reditoolsTable = pd.read_csv('REDItools2-Multiple_Table.csv').to_dict(orient='list')
# SPRINT
sprintTable = pd.read_csv('SPRINT-Multiple_Table.csv').to_dict(orient='list')
# RED-ML
redmlTable = pd.read_csv('REDML-Multiple_Table.csv').to_dict(orient='list')
# BCFtools
bcftoolsTable = pd.read_csv('BCFTools-Multiple_Table.csv').to_dict(orient='list')
# JACUSA
jacusaTable = pd.read_csv('JACUSA2_Table.csv').to_dict(orient='list')

# Plot within JACUSA2
create_plot_for_Jacusa  (prepare_for_jacusa(jacusaTable, aligners=aligners[0]), 'Jacusa2', aligners[0], 'JACUSA2')

############### Plot to compare 
tmp = {condition:{category:
        {align:[] for align in aligners[0]} 
        for category in categories} for condition in conditions}
AllData = {'BCFTools': bcftoolsTable, 
           'SPRINT':sprintTable, 
           'RED-ML': redmlTable, 
           'REDItools2':reditoolsTable,
            'JACUSA2':jacusaTable}
# Prepare for plotting
for tool in tools:
    data = AllData[tool]
    for condition in conditions:
        for align in aligners[0]: 
            if tool == 'SPRINT' and align == 'star':
                tmp[condition]['N_res'][align].append(0)
                tmp[condition]['N_db'][align].append(0)
                tmp[condition]['N_Alu'][align].append(0)
                continue
            # Idx for the data
            idxA = list(locate(data['Aligner'], lambda x: x == align.upper()))
            idxB = list(locate(data['Sample Condition'], lambda x: x == condition))
            idx = [i for i in idxA if i in idxB]
            # Save data to plot
            tmp[condition]['N_res'][align].append([x[1] for x in enumerate(data['# RES']) if x[0] in idx][0])
            tmp[condition]['N_db'][align].append([x[1] for x in enumerate(data['% RES in REDIportal']) if x[0] in idx][0])
            tmp[condition]['N_Alu'][align].append([x[1] for x in enumerate(data['% RES in Alu']) if x[0] in idx][0])

label_axisX = 'Tools' 
values_axisX = tools

# Plotting
def plot_case (data, posx, posy, c, multiplier=0, cont =0):
    for atribute, measurement in data.items():
        offset = width * multiplier - 0.2 # Location of bars
        rects = ax[posx,posy].bar(x + offset , measurement, width, label= aligners[0][cont], zorder=3, color= colors[cont])
        ax[posx,posy].set_xticks(x, subplot_xlabel, fontsize =15)
        ax[posx, posy].set_title(subplot_title[c], fontsize=20)
        ax[posx, posy].set_ylabel(subplot_ylabel[c], fontsize=17)
        ax[posx, posy].grid(zorder=0, axis='y')
        multiplier += 1 
        cont +=1  
    return

# Module 1 - Style format
golden_mean = (np.sqrt(5)-1.0)/2.0      # Aesthetic ratio
fig_width = 4                           # width in inches
fig_height = fig_width*golden_mean      # height in inches (4.5 for two vertical figures?)
fig_size = [fig_width,fig_height]
# Define format parameters
paramsscreen = {'backend': 'ps',
          'font.family':'serif',
          'axes.labelsize': 15,
           'legend.fontsize': 15,
           'xtick.labelsize': 12,
           'ytick.labelsize': 12,
           'figure.figsize': np.array(fig_size)}
paramsprint = {'backend': 'ps',
          'font.family':'serif',
          'axes.labelsize': 10,
            'legend.fontsize': 15,
           'xtick.labelsize': 8,
           'ytick.labelsize': 8,
           'figure.figsize': fig_size}
plt.rcParams.update(paramsscreen)

# Create the figure
fig, ax = plt.subplots(3,2,  figsize=(21, 18), dpi=300, sharex=True, sharey='row')
fig.text(0.3, 0.94, 'WT condition', ha='center', va='center', fontsize=30, fontweight='bold') # xlabel  
fig.text(0.7, 0.94, 'ADAR1KO condition', ha='center', va='center', fontsize=30, fontweight='bold') # xlabel 
fig.text(0.5, 0.08, label_axisX, ha='center', va='center', fontsize=20) # xlabel
# Module 2 - Plotting 
### Local Variables
subplot_title = ['# RES', 
                 '% RES in REDIportal', 
                 '% RES in Alu']
subplot_ylabel = ['Total counts', 
                  'Percentage (%)', 
                  'Percentage (%)']
colors = ['tomato',  # BWA
          'gold',  # HISAT2
          'dodgerblue'] # STAR
subplot_xlabel, xlab =  values_axisX, label_axisX  
x, width, c = np.arange(5), 0.2, 0

# Plot
present = 0
for condition in tmp:
    posx, c =0, 0
    posy = 0
    if condition == 'ADAR1KO':
        posy=1
    for case in tmp[condition]:
        plot_case (tmp[condition][case], posx, posy, c=c)
        posx +=1 
        c +=1
        present +=1

        
# Module 3 - Format plot 
axs = ax.flat
for n, ax in enumerate(axs): # Enumerate subplots as A, B, C ...
    ax.text(-0.07, 1.05, string.ascii_uppercase[n], transform=ax.transAxes, 
            size=20, weight='bold')
    ax.yaxis.set_tick_params(labelleft=True)
    labely=-0.2
    ax.spines['left'].set_position(('outward',10))
    ax.spines['bottom'].set_position(('outward',10))
    ax.spines['top'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

# Legend
lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
limit = 3
fig.legend(lines[0:limit], [label.upper() for label in labels[0:limit]], loc='lower center', ncol=3, title = "Aligner", 
               bbox_to_anchor=[0.5, 0.02], bbox_transform=fig.transFigure, title_fontsize='20')
plt.savefig('MultipleCompare.png', bbox_inches='tight')

        