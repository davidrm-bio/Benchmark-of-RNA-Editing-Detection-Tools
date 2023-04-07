### Analysis with duplicates plots ###
# Modules
import matplotlib.pyplot as plt
import ujson as json
import os
import numpy as np
import pandas as pd
import string

# Functions
def load_csv_multiple(file):
    cols = ['Condition', 'Aligner', 'Filter', 'RES', 'FDR', 'REDIportal', 'Alu']
    data = pd.read_csv(file, delimiter=',', header=None, names=cols)
    return data

def process_data_multiple (data, redml=False, sprint=False, bcftool=False):
    # Local Variables
    conditions = ['WT', 'ADAR1KO']
    aligners = ['BWA', 'STAR', 'HISAT2']
    supports = [2,4,6,8,10]
    categories = ['RES', 'REDIportal', 'Alu', 'FDR']
    
    if sprint: 
        aligners = ['BWA', 'HISAT2']
    if bcftool: 
        supports = [0, 0.1]
    if redml: 
        supports = [0.5,0.6,0.7,0.8,0.9]
        
    # Dictionary 
    data_plot = {condition: {category: {align: [] for align in aligners} 
                        for category in categories} for condition in conditions}

    # Separate data
    ycondition = data.Condition.to_list()
    yaligner =  data.Aligner.to_list()
    xfilter = data.Filter.to_list()
    yRES = data.RES.to_list()
    yFDR = data.FDR.to_list()
    yREDIportal = data.REDIportal.to_list()
    yAlu = data.Alu.to_list()
    
    # Save data
    for align in aligners:
        for condition in conditions:
            for support in supports:
                for idx in range (len (ycondition)):
                    if yaligner[idx] == align and ycondition[idx] == condition and xfilter[idx] == support:
                        data_plot[condition]['RES'][align].append (yRES[idx])
                        data_plot[condition]['REDIportal'][align].append (yREDIportal[idx])
                        data_plot[condition]['FDR'][align].append (yFDR[idx])
                        data_plot[condition]['Alu'][align].append (yAlu[idx])
                    
    return data_plot


def plot_within_tool_multiple(data, filename,  sprint=False, bcftool=False):
    # Plot
    fig, ax = plt.subplots(4,2, figsize=(21, 18), dpi=300, sharex=True, sharey='row')
    
    # Global Variables 
    subplot_title = ['# RES', '% RES in REDIportal', '% RES in Alu',  'FDR']
    subplot_ylabel = ['Total counts', 'Percentage (%)', 'Percentage (%)', 'Value']
    colors = ['#DCDCDC', '#BEBEBE', '#808080']
    xlab = 'Number of supporting reads'
    subplot_xlabel = [2,4,6,8,10]
    x, width, c = np.arange(5), 0.2, 0

    if bcftool: 
        xlab = 'Minor Allele Frequency (MAF)'
        subplot_xlabel = [0, 0.1]
        x = np.arange(2)
    
    fig.text(0.5, 0.093, xlab, ha='center', va='center', fontsize=20) # xlabel
    fig.text(0.3, 0.93, 'WT condition', ha='center', va='center', fontsize=30, fontweight='bold') 
    fig.text(0.7, 0.93, 'ADAR1KO condition', ha='center', va='center', fontsize=30, fontweight='bold') 
    
    # Plotting
    for condition in data:
        posx, posx2, c =0, 0, 0
        for case in data[condition]:
            cont = 0
            def plot_case (data, posx, posy, multiplier=0, cont =0, c=c):
                for atribute, measurement in data.items():
                    offset = width * multiplier - 0.2 # Location of bars
                    rects = ax[posx,posy].bar(x + offset , measurement, width, label= aligners[cont], zorder=3, color= colors[cont] )
                    ax[posx,posy].set_xticks(x, subplot_xlabel, fontsize =15)
                    ax[posx, posy].set_title(subplot_title[c], fontsize=20)
                    ax[posx, posy].set_ylabel(subplot_ylabel[c], fontsize=17)
                    ax[posx, posy].grid(zorder=0, axis='y')
                    multiplier += 1 
                    cont +=1  
                return
            if  condition == 'WT':
                plot_case (data['WT'][case], posx, 0)
            else:
                plot_case (data['ADAR1KO'][case], posx, 1)
            posx +=1 
            c +=1

    # General setting post-plot  
    axs = ax.flat
    for n, ax in enumerate(axs): # Enumerate subplots as A, B, C ...
        ax.text(-0.07, 1.05, string.ascii_uppercase[n], transform=ax.transAxes, 
                size=20, weight='bold')
        ax.yaxis.set_tick_params(labelleft=True)

    lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
    lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
    limit = 3
    if sprint: 
        limit = 2
        labels = ['BWA', 'HISAT2']
    
    fig.legend(lines[0:limit], labels[0:limit], loc='lower center', ncol=3, prop={'size': 13}, title = "Aligner", 
               bbox_to_anchor=[0.5, 0.03], bbox_transform=fig.transFigure)
    
    # Save fig
    plt.savefig(filename  + '.png', bbox_inches='tight')
    return



def plot_within_redml_multiple(data, filename):
    # Plot
    fig, ax = plt.subplots(3,2,  figsize=(21, 18), dpi=300, sharex=True, sharey='row')
    
    fig.text(0.5, 0.05, 'Detection Threshold', ha='center', va='center', fontsize=17) # xlabel
    fig.text(0.3, 0.93, 'WT condition', ha='center', va='center', fontsize=30, fontweight='bold') 
    fig.text(0.7, 0.93, 'ADAR1KO condition', ha='center', va='center', fontsize=30, fontweight='bold') 

    # Global Variables 
    subplot_title = ['# RES', '% RES in REDIportal', '% RES in Alu']
    subplot_ylabel = ['Total counts', 'Percentage (%)', 'Percentage (%)', 'Value']
    colors = ['#DCDCDC', '#BEBEBE', '#808080']
    subplot_xlabel = [0.5,0.6,0.7,0.8,0.9]
    x, width = np.arange(5), 0.2

    for condition in data:
        posx, c = 0, 0
        for case in data[condition]:
            if 'FDR'  not in case:
                cont = 0
                def plot_case (data, posx, posy, multiplier=0, cont =0, c=c):
                    for atribute, measurement in data.items():
                        offset = width * multiplier - 0.2 # Location of bars
                        rects = ax[posx,posy].bar(x + offset , measurement, width, label= aligners[cont], zorder=3, color= colors[cont] )
                        ax[posx,posy].set_xticks(x, subplot_xlabel, fontsize =15)
                        ax[posx, posy].set_title(subplot_title[c], fontsize=20)
                        ax[posx, posy].set_ylabel(subplot_ylabel[c], fontsize=17)
                        ax[posx, posy].grid(zorder=0, axis='y')
                        multiplier += 1 
                        cont +=1  
                    return
                if  condition == 'WT':
                    plot_case (data['WT'][case], posx, 0)
                else:
                    plot_case (data['ADAR1KO'][case], posx, 1)

                posx +=1 
                c +=1

    # General setting post-plot  
    axs = ax.flat
    for n, ax in enumerate(axs): # Enumerate subplots as A, B, C ...
        ax.text(-0.07, 1.05, string.ascii_uppercase[n], transform=ax.transAxes, 
                size=20, weight='bold')
        ax.yaxis.set_tick_params(labelleft=True)

    lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
    lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
    fig.legend(lines[0:3], labels[0:3], loc='upper center', ncol=3, prop={'size': 13}, title = "Aligner")

    # Save fig
    plt.savefig(filename  + '.png', bbox_inches='tight')
    return    

# Global variables
aligners = ['BWA', 'STAR', 'HISAT2'] # ['bwa', 'hisat2']
conditions = ['WT', 'ADAR1KO']
tools = ['BCFtools', 'REDML', 'SPRINT', 'REDItools2', 'JACUSA2']


# Load data
BCFtools = load_csv_multiple('BCFtools_multiple_Excel.csv')
REDML = load_csv_multiple('REDML_multiple_Excel.csv')
SPRINT = load_csv_multiple('SPRINT_multiple_Excel.csv')
REDItools2 = load_csv_multiple('REDItools2_multiple_Excel.csv')
Jacusa2 = load_csv_multiple ('Jacusa2_Excel.csv')

data_plot = {tool: {} for tool in tools}

bcftools_plot = process_data_multiple(BCFtools, bcftool=True)
redml_plot = process_data_multiple(REDML, redml=True)
sprint_plot = process_data_multiple(SPRINT, sprint=True)
reditools2_plot = process_data_multiple(REDItools2)
jacusa2_plot = process_data_multiple(Jacusa2)

data_plot['BCFtools'] = bcftools_plot
data_plot['REDML'] = redml_plot
data_plot['SPRINT'] = sprint_plot
data_plot['REDItools2'] = reditools2_plot
data_plot['JACUSA2'] = jacusa2_plot


plot_within_tool_multiple (data_plot['REDItools2'], 'REDItools2_multiple')
plot_within_tool_multiple (data_plot['BCFtools'], 'BCFtools_individual', bcftool=True)
plot_within_redml_multiple (data_plot['REDML'], 'REDML_multiple')
plot_within_tool_multiple (data_plot['SPRINT'], 'SPRINT_multiple',  sprint=True)
plot_within_tool_multiple (data_plot['JACUSA2'], 'JACUSA2_multiple')

# DATA to comapre multiple
tools = ['BCFtools', 'REDML', 'SPRINT', 'REDItools2', 'JACUSA2']
conditions = ['WT', 'ADAR1KO']
aligners = ['BWA', 'STAR', 'HISAT2']
categories = ['RES', 'REDIportal', 'Alu', 'FDR']

for condition in conditions:
    for category in categories:
        sprint_plot[condition][category]['STAR'] = [0]

data_compare = {condition: {category: {align: [] for align in aligners} for category in categories} for condition in conditions}

for idx, tool in enumerate ([bcftools_plot, redml_plot, sprint_plot, reditools2_plot, jacusa2_plot]):
    for condition in conditions:
        for category in categories:
            for align in aligners:
                data_compare[condition][category][align].append(tool[condition][category][align][0])
                
data = data_compare

# Plot
fig, ax = plt.subplots(4,2, figsize=(21, 18), dpi=300, sharex=True, sharey='row')

# Global Variables 
subplot_title = ['# RES', '% RES in REDIportal', '% RES in Alu',  'FDR']
subplot_ylabel = ['Total counts', 'Percentage (%)', 'Percentage (%)', 'Value']
subplot_xlabel = ['BCFtools', 'REDML', 'SPRINT', 'REDItools2', 'JACUSA2']
colors = ['#DCDCDC', '#BEBEBE', '#808080']
posx, posy = 0,0
x, width, c = np.arange(5), 0.2, 0

fig.text(0.5, 0.093, 'Benchmarked tool', ha='center', va='center', fontsize=20) # xlabel
fig.text(0.3, 0.93, 'WT condition', ha='center', va='center', fontsize=30, fontweight='bold') 
fig.text(0.7, 0.93, 'ADAR1KO condition', ha='center', va='center', fontsize=30, fontweight='bold') 

# Plotting
for condition in data:
    posx, posx2, c =0, 0, 0
    for case in data[condition]:
        cont = 0
        def plot_case (data, posx, posy, multiplier=0, cont =0, c=c):
            for atribute, measurement in data.items():
                offset = width * multiplier - 0.2 # Location of bars
                rects = ax[posx,posy].bar(x + offset , measurement, width, label= aligners[cont], zorder=3, color= colors[cont] )
                ax[posx,posy].set_xticks(x, subplot_xlabel, fontsize =15)
                ax[posx, posy].set_title(subplot_title[c], fontsize=20)
                ax[posx, posy].set_ylabel(subplot_ylabel[c], fontsize=17)
                ax[posx, posy].grid(zorder=0, axis='y')
                multiplier += 1 
                cont +=1  
            return
        if  condition == 'WT':
            plot_case (data['WT'][case], posx, 0)
        else:
            plot_case (data['ADAR1KO'][case], posx, 1)
        posx +=1 
        c +=1

# General setting post-plot  
axs = ax.flat
for n, ax in enumerate(axs): # Enumerate subplots as A, B, C ...
    ax.text(-0.07, 1.05, string.ascii_uppercase[n], transform=ax.transAxes, 
            size=20, weight='bold')
    ax.yaxis.set_tick_params(labelleft=True)

lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
limit = 3
fig.legend(lines[0:limit], labels[0:limit], loc='lower center', ncol=3, prop={'size': 13}, title = "Aligner", 
           bbox_to_anchor=[0.5, 0.03], bbox_transform=fig.transFigure)
# Save fig
plt.savefig('Compare_multiple'  + '.png', bbox_inches='tight')

