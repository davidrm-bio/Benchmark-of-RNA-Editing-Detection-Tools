# Modules
import os
import pandas as pd
import numpy as np
import ujson as json
import warnings
from matplotlib import rc                 #To alter the defaults of plots       
from IPython import display   
import string
import matplotlib.pyplot as plt
from more_itertools import locate
import matplotlib as mpl
warnings.filterwarnings("ignore")


# Functions


def load_file (filename):
    '''Load tab-delimited files into a DataFrame'''
    df = pd.read_csv(os.path.join(filename), delimiter="\t") 
    return df

def load_database (path, file): 
    '''Load JSON file into dictionary '''
    with open (os.path.join(path, file), 'r') as file:
        db = json.load (file)
    return db


def count_Alu (data, alu):
    '''
    providing a list of RES calculate
    how many are in Alu regions
    '''
    cont =0
    for IDs in data:
        chrom = IDs.split('|')[0]
        pos = IDs.split('|')[1]
        if chrom in alu:
            if pos in alu[chrom]:
                cont += 1
    return cont


def save2JSON (path, filename, data):
    with open (os.path.join (path, filename), 'w') as out:
        json.dump(data, out)
    return


def stats (filename, tool, support, rediportal):
    '''
    Compute: Number  of RES; Number of RES in REDIportal; Number of RES in Alu & SNPs
    '''
    with open (filename, 'r')  as Read:
        # Get all RES
        list_of_res = []
        for lines in Read:
            if '#' in lines:
                SNPs = lines.strip().split(' ')[4]
            else:
                line = lines.strip().split('\t')
                # Get Res depending on the  tool
                if tool == 'REDItools2':
                    if float(line[4]) >= support:
                        list_of_res.append(line[0] + '|' + line[1] + '|' + line[7]) # REDItools2
                elif tool == 'REDML':
                    if  float(line[-1]) >= support:
                        list_of_res.append(line[0] + '|' + line[1] + '|' + line[3] + line[5]) # REDML
                elif tool == 'SPRINT':
                    if  float(line[4]) >= support:
                        list_of_res.append(line[0] + '|' + line[2] + '|' + line[3]) # REDML
                elif tool == 'BCFTools':
                    list_of_res.append(line[0] + '|' + line[2] + '|' + line[5] + line[6]) # BCFtools
                else:
                    print ('Error, not an option')
    # Check with REDIportal
    res_in_db = [res for res in list_of_res if res in rediportal]
    return len(list_of_res), len(res_in_db), int(SNPs), list_of_res



def save_to_csv (filename, data, aligners, conditions, samples, thresholds, thresholdname='Number of Supporting Reads'):
    '''
    Export to CSV the data that is plotted
    '''
    line0 = ['Sample Condition', 
                'Aligner', 
                thresholdname, 
                'Average # RES', 
                'RES SEM', 
                'Average % RES in REDIportal', 
                'Average % RES in Alu',
                'DB SEM',
                'Alu SEM']
    # Save data
    with open (filename, 'w') as out:
        out.write (','.join(line0) + '\n')
        for align in aligners:
            for condition in conditions:
                for cutoff in thresholds:
                    clones_res = [data[align][condition][clone][cutoff]['N_RES'] for clone in samples]
                    Averg_RES = np.mean (clones_res)
                    RES_SEM = sem (clones_res)
                    Averg_db = np.mean ([data[align][condition][clone][cutoff]['N_REDIportal']/data[align][condition][clone][cutoff]['N_RES']*100 for clone in samples])
                    Averg_alu = np.mean ([data[align][condition][clone][cutoff]['N_Alu']/data[align][condition][clone][cutoff]['N_RES']*100 for clone in samples])
                    # Extra needed for plot
                    sem_db = sem ([data[align][condition][clone][cutoff]['N_REDIportal']/data[align][condition][clone][cutoff]['N_RES']*100 for clone in samples])
                    sem_alu = sem([data[align][condition][clone][cutoff]['N_Alu']/data[align][condition][clone][cutoff]['N_RES']*100 for clone in samples])
                    newline = [align.upper(), 
                            condition, 
                            str(cutoff), 
                            str(int (round (Averg_RES,0))),
                            str(round (RES_SEM,2)),
                            str(round (Averg_db,2)),
                            str(round (Averg_alu,2)),
                            str(round (sem_db,2)),
                            str(round (sem_alu,2))]
                    out.write(','.join(newline) + '\n')
    return


def save_to_csv_multiple (filename, data, aligners, conditions, thresholds, thresholdname='Number of Supporting Reads'):
    '''
    Export to CSV the data that is plotted - Multiple Analysis
    '''
    line0 = ['Aligner', 
                'Sample Condition', 
                thresholdname, 
                '# RES', 
                '% RES in REDIportal', 
                '% RES in Alu']
    # Save data
    with open (filename, 'w') as out:
        out.write (','.join(line0) + '\n')
        for align in aligners:
            for condition in conditions:
                for cutoff in thresholds:
                    clones_res = data[align][condition][cutoff]['N_res']
                    try: 
                        Averg_db = data[align][condition][cutoff]['N_db']/data[align][condition][cutoff]['N_res']*100
                    except ZeroDivisionError:
                        Averg_db = 0 # SPRINT for multiple has cases with 0 in the ADAR1KO condition
                    try: 
                        Averg_alu = data[align][condition][cutoff]['N_Alu']/data[align][condition][cutoff]['N_res']*100
                    except ZeroDivisionError:
                        Averg_alu = 0  # SPRINT for multiple has cases with 0 in the ADAR1KO condition
                        
                    # Extra needed for plot
                    newline = [align.upper(), 
                            condition, 
                            str(cutoff), 
                            str(int (clones_res)),
                            str(round (Averg_db,2)),
                            str(round (Averg_alu,2))]
                    out.write(','.join(newline) + '\n')
    return

def generate_snp_table(data, aligners, conditions, samples, toolname, cutoff='2', firstline=False, filename='IndividualAnalysis-SNPs-Table.csv'):
    '''
    Generate CSV that summaries the number of SNPs detected by each tool
    '''
    line0 = ['Tool',
             'Sample Condition', 
            'Aligner', 
            'Average # SNPs', 
            'SNPs SEM']
    with open (filename, 'a') as out:
        if firstline:
            out.write(','.join(line0) + '\n')
        for align in aligners:
            for condition in conditions:
                for clone in samples:
                    if toolname == 'JACUSA2':
                        Averg_snps = data[align][condition][cutoff]['SNPs']
                        SEM_snps = 0
                    else:
                        clones_snps = [data[align][condition][clone][cutoff]['SNPs'] for clone in samples]
                        Averg_snps = np.mean (clones_snps)
                        SEM_snps = sem (clones_snps)
                newline = [toolname, 
                           condition, 
                            align.upper(), 
                            str(int (round (Averg_snps, 0))),
                            str(round (SEM_snps,2))]
                out.write(','.join(newline) + '\n')
    return


def prepare_for_plotting (data, aligners, conditions = ['WT', 'ADAR1KO']):
    # Local Variables
    categories = ['N_res', 'SEM_res', 'N_db', 'SEM_db', 'N_Alu', 'SEM_Alu']
    tmp = {condition: {category:
           {align:[] for align in aligners} 
           for category in categories} for condition in conditions}
    # Analysis
    for condition in conditions:
        for align in aligners:
            # Idx for the data
            idxA = list(locate(data['Sample Condition'], lambda x: x == align.upper()))
            idxB = list(locate(data['Aligner'], lambda x: x == condition))
            idx = [i for i in idxA if i in idxB]
            # Save data to plot
            tmp[condition]['N_res'][align] = [x[1] for x in enumerate(data['Average # RES']) if x[0] in idx]
            tmp[condition]['SEM_res'][align] = [x[1] for x in enumerate(data['RES SEM']) if x[0] in idx]

            tmp[condition]['N_db'][align] = [x[1] for x in enumerate(data['Average % RES in REDIportal']) if x[0] in idx]
            tmp[condition]['SEM_db'][align] = [x[1] for x in enumerate(data['DB SEM']) if x[0] in idx]

            tmp[condition]['N_Alu'][align] = [x[1] for x in enumerate(data['Average % RES in Alu']) if x[0] in idx]
            tmp[condition]['SEM_Alu'][align] = [x[1] for x in enumerate(data['Alu SEM']) if x[0] in idx]
    return tmp


def create_plot (tmp, filename, aligners, title, label_axisX = 'Number of supporting reads' , values_axisX = [2,4,6,8,10],  sprint = False, bcftools=False):
    
    # Internal functions
    def plot_case (data, posx, posy, c, multiplier=0, cont =0):
        for atribute, measurement in data.items():
            offset = width * multiplier - 0.2 # Location of bars
            rects = ax[posx,posy].bar(x + offset , measurement, width, label= aligners[cont], zorder=3, color= colors[cont])
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
    
    fig.text(0.5, 0.98, 'Statistics on ' + title, ha='center', va='center', fontsize=35, fontweight='bold') # Title
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
    x, width, c = np.arange(5), 0.2, 0
    if bcftools:
        x = np.arange(2)
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
    if sprint:
        limit = 2
    fig.legend(lines[0:limit], [label.upper() for label in labels[0:limit]], loc='lower center', ncol=3, title = "Aligner", 
                   bbox_to_anchor=[0.5, 0.02], bbox_transform=fig.transFigure, title_fontsize='20')
    plt.savefig(filename  + '.png', bbox_inches='tight')
    return


def Generate_Json_Multiple(path, filename, data, rediportal, alu, aligners, thresholds):
    '''
    Re-Analysis to keep RES present in the three replicates
    '''
    # Local Variables
    categories = ['N_res', 'N_db', 'N_Alu', 'List_of_Res']
    conditions = ['WT', 'ADAR1KO']
    export_data = {align: 
                   {condition:
                    {cutoff:
                     {category:0 for category in categories}
                    for cutoff in thresholds} 
                    for condition in conditions} 
                  for align in aligners}
    # Generate Data
    for align in aligners:          
        for condition in conditions:
            for cutoff in thresholds:
                clone1 = data[align][condition]['clone1'][cutoff]['List_of_Res']
                clone2 = data[align][condition]['clone2'][cutoff]['List_of_Res']
                clone3 = data[align][condition]['clone3'][cutoff]['List_of_Res']
                # Res present in the 3 replicates
                common = set(clone1) & set(clone2) & set(clone3)
                # Support in REDIportal
                N_db = [res for res in common if res in rediportal]
                # Count in Alu regions
                N_alu = count_Alu(common, alu)
                # Save data to export
                export_data[align][condition][cutoff]['N_res'] = len (common)
                export_data[align][condition][cutoff]['N_db'] = len(N_db)
                export_data[align][condition][cutoff]['N_Alu'] = N_alu
                export_data[align][condition][cutoff]['List_of_Res'] = list(common)
    # Save data
    save2JSON (path, filename, export_data)
    return


def prepare_for_jacusa(jacusaTable, aligners, conditions = ['WT', 'ADAR1KO'] ):
    # Local variables
    categories = ['N_res', 'N_db', 'N_Alu']
    tmp = {condition: {category:
           {align:[] for align in aligners} 
           for category in categories} for condition in conditions}
    # Analysis
    for condition in conditions:
        for align in aligners:
            # Idx for the data
            idxA = list(locate(data['Aligner'], lambda x: x == align.upper()))
            idxB = list(locate(data['Sample Condition'], lambda x: x == condition))
            idx = [i for i in idxA if i in idxB]
            # Save data to plot
            tmp[condition]['N_res'][align] = [x[1] for x in enumerate(data['# RES']) if x[0] in idx]
            tmp[condition]['N_db'][align] = [x[1] for x in enumerate(data['% RES in REDIportal']) if x[0] in idx]
            tmp[condition]['N_Alu'][align] = [x[1] for x in enumerate(data['% RES in Alu']) if x[0] in idx]
    return tmp

def create_plot_for_Jacusa (tmp, filename, aligners, title, label_axisX = 'Number of supporting reads' , values_axisX = [2,4,6,8,10]):
    # Internal functions
    def plot_case (data, posx, posy, c, multiplier=0, cont =0):
        for atribute, measurement in data.items():
            offset = width * multiplier - 0.2 # Location of bars
            rects = ax[posx,posy].bar(x + offset , measurement, width, label= aligners[cont], zorder=3, color= colors[cont])
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
    
    fig.text(0.5, 0.98, 'Statistics on ' + title, ha='center', va='center', fontsize=35, fontweight='bold') # Title
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
        ax.yaxis.set_major_formatter(mpl.ticker.StrMethodFormatter('{x:,.0f}'))

    # Legend
    lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
    lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
    limit = 3
    fig.legend(lines[0:limit], [label.upper() for label in labels[0:limit]], loc='lower center', ncol=3, title = "Aligner", 
                   bbox_to_anchor=[0.5, 0.02], bbox_transform=fig.transFigure, title_fontsize='20')
    plt.savefig(filename  + '.png', bbox_inches='tight')
    return
