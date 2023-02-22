#!/usr/bin/env python

import glob
files = glob.glob('*.bam')
files_processed = [name.replace('.bam', '') for name in files]
with open('samples.txt', 'w') as f:
	f.write('\n'.join(map(str, files_processed)))
    