# Modules
import os
import ujson as json
import pandas as pd
import warnings

warnings.filterwarnings("ignore")

# Define Global variables
path = './'
files = os.listdir(path)
vcf_files = [file for file in files if 'vcf' in file]
cols = [0,1,3,4] # [chr, pos, ref, alt]


for file in vcf_files:
    # Define Local variables
    filename = file.replace('vcf', 'json')
    
    # Load data
    vcf = pd.read_csv(os.path.join (path, file), sep="\t", comment='#', header=None)

    # Select columns to generate ID
    vcf_filt = vcf[cols]

    # Filter to retain only SNV
    vcf_filt['Ref_len'] = vcf_filt[3].str.len()
    vcf_filt['Alt_len'] = vcf_filt[4].str.len()

    # Only SNVs
    vcf_filt_snp = vcf_filt[(vcf_filt['Ref_len'] == 1) & ( vcf_filt['Alt_len']==1)]

    # Create ID Chr1|1989723|AG
    vcf_filt_snp['ID'] = vcf_filt_snp[0].astype(str) + '|' + vcf_filt_snp[1].astype(str) + '|' + vcf_filt_snp[3] + vcf_filt_snp[4]

    # Convert from Pandas > JSON
    vcf_dict = dict(zip(vcf_filt_snp.ID, vcf_filt_snp.ID))
    
    with open(filename, 'w') as fp:
        json.dump(vcf_dict, fp)