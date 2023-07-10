# Modules
import ujson as json

def bed2json (filename):
    '''
    Convert BED file to json format. 
    Takes the last column that contains the RES Ids
    '''
    newfilename = filename.replace('.bed', '.json')
    with open(newfilename, 'w') as out:
        with open(filename, 'r') as Read:
            # Get RES Ids
            ids = []
            for lines in Read:
                line = lines.strip().split('\t')
                ids.append(line[-1])
            # Convert to dict
            data  = dict(zip(ids,ids))
            # Dict to Json
            json.dump(data, out)
    return

            
bed2json('HEK293T_hg19_clean.bed')
bed2json('HEK293T_hg38_clean.bed')
    
                
                