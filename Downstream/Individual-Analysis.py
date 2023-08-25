# Modules 
import ujson as json
from MainFunctions import  *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats as sem
from more_itertools import locate
import matplotlib as mpl

# Global variables
path = '..'
aligners = [['bwa', 'hisat2', 'star'],  # For the rest
            ['bwa', 'hisat2']] # For SPRINT
samples = ['clone1', 'clone2', 'clone3']
conditions = ['WT', 'ADAR1KO']
#categories = ['N_RES', 'N_REDIportal', 'N_Alu', 'SNPs', 'List_of_Res']
thresholds = [['2', '4', '6', '8', '10'],  # For the rest
             ['0.5','0.6', '0.7', '0.8', '0.9'], # For REDML
             ['0', '0.1']] # For BCFtools


# Load Data
reditools = load_database(path, 'Data_REDItool2.json')
redml = load_database(path, 'Data_REDML.json')
sprint = load_database(path, 'Data_SPRINT.json')
bcftools = load_database(path, 'Data_BCFTools.json')
# Generate Tables
save_to_csv('REDItools2_Table.csv', reditools, aligners[0], conditions, samples, thresholds[0])
save_to_csv('SPRINT_Table.csv', sprint, aligners[1], conditions, samples, thresholds[0])
save_to_csv('REDML_Table.csv', redml, aligners[0], conditions, samples, thresholds[1], thresholdname='Detection Threshold')
save_to_csv('BCFTools_Table.csv', bcftools, aligners[0], conditions, samples, thresholds[2], thresholdname = 'Minor Allele Frequency')
# Generate SNP Table
generate_snp_table(bcftools, aligners[0], conditions, samples, cutoff='0', firstline=True, toolname ='BCFTools')
generate_snp_table(redml, aligners[0], conditions, samples, cutoff='0.5', toolname ='RED-ML')
generate_snp_table(sprint, aligners[1], conditions, samples, cutoff='2', toolname ='SPRINT')
generate_snp_table(reditools, aligners[0], conditions, samples, cutoff='2', toolname ='REDItools2')
generate_snp_table(jacusa, aligners[0], conditions, samples, cutoff='2', toolname ='JACUSA2')

# Plots within tools

# REDItools2
reditoolsTable = pd.read_csv('REDItools2_Table.csv').to_dict(orient='list')
create_plot(prepare_for_plotting(reditoolsTable, aligners[0]), 'REDItools2_Individual', aligners=aligners[0], title='REDItools2')
# SPRINT
sprintTable = pd.read_csv('SPRINT_Table.csv').to_dict(orient='list')
create_plot(prepare_for_plotting(sprintTable, aligners[1]), 'SPRINT_Individual', sprint=True, aligners=aligners[1], title='SPRINT')
# RED-ML
redmlTable = pd.read_csv('REDML_Table.csv').to_dict(orient='list')
create_plot(prepare_for_plotting(redmlTable, aligners[0]), 'REDML_Individual', aligners=aligners[0], label_axisX = 'Detection Threshold' , values_axisX = [0.5,0.6,0.7,0.8,0.9], title='RED-ML')
# BCFtools
bcftoolsTable = pd.read_csv('BCFTools_Table.csv').to_dict(orient='list')
create_plot(prepare_for_plotting(bcftoolsTable, aligners[0]), 'BCFTools_Individual', aligners=aligners[0], label_axisX = 'Minor Allele Frequency' , values_axisX = [0,0.1], bcftools=True, title='BCFTools')


# Plot to compare - #################################  The script below was run on jupyter notebook - Using all modules loaded here
categories = ['N_res', 'SEM_res', 'N_db', 'SEM_db', 'N_Alu', 'SEM_Alu']
tools = ['BCFTools', 'RED-ML', 'SPRINT', 'REDItools2']

## Variables
tmp = {condition:{category:
        {align:[] for align in aligners[0]} 
        for category in categories} for condition in conditions}
AllData = {'BCFTools': bcftoolsTable, 
           'SPRINT':sprintTable, 
           'RED-ML': redmlTable, 
           'REDItools2':reditoolsTable}

label_axisX = 'Tools' 
values_axisX = tools

# Prepare for plotting
for tool in tools:
    data = AllData[tool]
    for condition in conditions:
        for align in aligners[0]: 
            if tool == 'SPRINT' and align == 'star':
                tmp[condition]['N_res'][align].append (0)
                tmp[condition]['SEM_res'][align].append(0)

                tmp[condition]['N_db'][align].append (0)
                tmp[condition]['SEM_db'][align].append(0)

                tmp[condition]['N_Alu'][align].append(0)
                tmp[condition]['SEM_Alu'][align].append(0)
                continue
            # Idx for the data
            idxA = list(locate(data['Sample Condition'], lambda x: x == align.upper()))
            idxB = list(locate(data['Aligner'], lambda x: x == condition))
            idx = [i for i in idxA if i in idxB]
            # Save data to plot
            tmp[condition]['N_res'][align].append ([x[1] for x in enumerate(data['Average # RES']) if x[0] in idx][0])
            tmp[condition]['SEM_res'][align].append([x[1] for x in enumerate(data['RES SEM']) if x[0] in idx][0])

            tmp[condition]['N_db'][align].append ([x[1] for x in enumerate(data['Average % RES in REDIportal']) if x[0] in idx][0])
            tmp[condition]['SEM_db'][align].append([x[1] for x in enumerate(data['DB SEM']) if x[0] in idx][0])

            tmp[condition]['N_Alu'][align].append([x[1] for x in enumerate(data['Average % RES in Alu']) if x[0] in idx][0])
            tmp[condition]['SEM_Alu'][align].append([x[1] for x in enumerate(data['Alu SEM']) if x[0] in idx][0])

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
def plot_errorbar (datax, datay, posx, posy,c,  multiplier2 =0):
    for atribute, measurement in datay.items():
        offset = width * multiplier2 - 0.2 # Location of bars
        ax[posx,posy].errorbar(x + offset, datax[atribute] , measurement, fmt ='o', elinewidth = 3, capsize=7, zorder=5, color ='black')
        multiplier2 += 1 
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
subplot_title = ['Average # RES', 
                 'Average % RES in REDIportal', 
                 'Average % RES in Alu']
subplot_ylabel = ['Total counts', 
                  'Percentage (%)', 
                  'Percentage (%)']
colors = ['tomato',  # BWA
          'gold',  # HISAT2
          'dodgerblue'] # STAR
subplot_xlabel, xlab =  values_axisX, label_axisX  
x, width, c = np.arange(4), 0.2, 0

# Plot
present = 0
for condition in tmp:
    posx, posx2, c =0, 0, 0
    posy = 0
    if condition == 'ADAR1KO':
        posy=1
    for case in tmp[condition]:
        if 'SEM' not in case:
            plot_case (tmp[condition][case], posx, posy, c=c)
            posx +=1 
            c +=1
        else:
            plot_errorbar(tmp[condition][case.replace('SEM', 'N')], tmp[condition][case], posx2, posy, c=c)
            posx2 += 1
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
    ax.yaxis.set_major_formatter(mpl.ticker.StrMethodFormatter('{x:,.0f}'))

# Legend
lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
limit = 3
fig.legend(lines[0:limit], [label.upper() for label in labels[0:limit]], loc='lower center', ncol=3, title = "Aligner", 
               bbox_to_anchor=[0.5, 0.02], bbox_transform=fig.transFigure, title_fontsize='20')
plt.savefig('IndividualCompare.png', bbox_inches='tight')
