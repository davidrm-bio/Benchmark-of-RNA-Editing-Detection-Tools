
conditions = ['WT', 'ADAR1KO']
aligners = ['BWA', 'HISAT2', 'STAR']
categories = ['N_SNPs', 'SEM_SNPs']
tools = ['BCFTools', 'RED-ML','SPRINT', 'REDItools2', 'JACUSA2']

#tmp = {condition: {category:
#       {align:{tool:[]for tool in tools} for align in aligners} 
#       for category in categories} for condition in conditions}

tmp = {condition: {category:
       {tool:[]for tool in tools}
       for category in categories} for condition in conditions}


# Plot number of SNPs reported by each tool
# Modules
import os
import pandas as pd
import numpy as np
import ujson as json
import warnings
import matplotlib.pyplot as plt
from matplotlib import rc                 #To alter the defaults of plots       
from IPython import display   
import string
from more_itertools import locate
import matplotlib as mpl
warnings.filterwarnings("ignore")

for condition in conditions:
    for align in aligners:
        for tool in tools:
            if tool == 'SPRINT' and align == 'STAR':
                tmp[condition]['N_SNPs'][tool].append(0)
                tmp[condition]['SEM_SNPs'][tool].append(0)
                continue
            # Idx for the data
            idxA = list(locate(data['Aligner'], lambda x: x == align.upper()))
            idxB = list(locate(data['Sample Condition'], lambda x: x == condition))
            idxC = list(locate(data['Tool'], lambda x: x == tool))
            idx = [i for i in idxA if i in idxB and i in idxC]
            # Save data to plot
            tmp[condition]['N_SNPs'][tool].append([x[1] for x in enumerate(data['Average # SNPs']) if x[0] in idx][0])
            tmp[condition]['SEM_SNPs'][tool].append([x[1] for x in enumerate(data['SNPs SEM']) if x[0] in idx][0])

# Variables
width = 1
posx, posx2, c =0, 0, 0
x = np.arange(0,20,7)
xlabels = ['BWA', 'HISAT2', 'STAR']
subplot_title = ['WT', 'ADAR1KO']
colors =['tomato', 'gold', 'palegreen', 'dodgerblue', 'darkorchid']
def plot_errorbar (datax, datay, posx,c,  multiplier2 =0):
    pos_report = []
    for atribute, measurement in datay.items():
        offset = width * multiplier2 - 0.2 # Location of bars
        ax[posx].errorbar(x + offset, datax[atribute] , measurement, fmt ='o', linewidth = 3, capsize=7, zorder=5, color ='black')
        multiplier2 += 1 
        pos_report.append(measurement)
    return pos_report

# Define format parameters

# Module 1 - Style format
golden_mean = (np.sqrt(5)-1.0)/2.0      # Aesthetic ratio
fig_width = 20                         # width in inches
fig_height = fig_width*golden_mean      # height in inches (4.5 for two vertical figures?)
fig_size = [fig_width,fig_height]
label_axisX = 'Tools'
values_axisX = tools


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

# Create plot
fig, ax = plt.subplots(1,2, dpi=300, sharex=True, sharey='row')
fig.text(0.5, 0.99, 'SNPs reported', ha='center', va='center', fontsize=35, fontweight='bold') # Title
fig.text(0.5, 0.06, 'Aligners', ha='center', va='center', fontsize=20) # xlabel

# Plotting

# Module 2 - Plotting 
for condition in tmp:
    multiplier=0
    cont = 0
    for case in tmp[condition]:
        if 'SEM' not in case:
            for atribute, measurement in tmp[condition][case].items():
                offset = width * multiplier - 0.2 # Location of bars
                rects = ax[posx].bar(x + offset , measurement, width, zorder=3, color= colors[cont], label = atribute)
                multiplier += 1 
                cont +=1 
            ax[posx].set_xticks(x + width, xlabels, fontsize =15)
            ax[posx].set_title(subplot_title[c], fontsize=20)
            ax[posx].set_ylabel('Total Count', fontsize=17)
            ax[posx].grid(zorder=0, axis='y')
            
        else:
            pos_report = plot_errorbar(tmp[condition][case.replace('SEM', 'N')], tmp[condition][case], posx2, c=c)

            i = 0
            j = 0
            for rect in ax[posx].patches:
                # Get X and Y placement of label from rect.
                y_value = rect.get_height() 
                x_value = rect.get_x() + rect.get_width() / 2
                # Number of points between bar and label. Change to your liking.
                space = pos_report[j][i]
                i += 1
                #space = 5
                # Vertical alignment for positive values
                va = 'top'
                # Use Y value as label and format number with one decimal place
                label = "{:,.0f}".format(y_value)
                # Create annotation
                ax[posx].annotate(
                    label,                      # Use `label` as label
                    (x_value, y_value+space),         # Place label at end of the bar
                    xytext=(0, 10),          # Vertically shift label by `space`
                    textcoords="offset points", # Interpret `xytext` as offset in points
                    ha='center',                # Horizontally center label
                    va=va)                      # Vertically align label differently for
                                                    # positive and negative values.
                if i == 3:
                    i = 0
                    j +=1
                
    posx += 1
    posx2 += 1
    c +=1

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
limit = 5
fig.legend(lines[0:limit], [label.upper() for label in labels[0:limit]], loc='lower center', ncol=3, title = "Tools", 
               bbox_to_anchor=[0.5, -0.07], bbox_transform=fig.transFigure, title_fontsize='20')
plt.savefig('SNPs Reported'  + '.png', bbox_inches='tight')

