#!/usr/bin/env python
# coding: utf-8

import glob
import sys
import os
import pandas as pd

def filter_pathway(cutoff, indir, outdir):
    '''
    Filter pathways whose completeness < cutoff%
    '''
    os.makedirs(outdir, exist_ok=True)
    for file in glob.glob(os.path.join(indir, '*.xlsx')):
        if '$' not in file:
            df = pd.read_excel(file)
            df = df[df['Pathway_completeness'] >= float(cutoff)] 
            df.to_excel(os.path.join(outdir, os.path.basename(file)), index=False)

if __name__ == '__main__':
    cutoff = sys.argv[1]
    mag_dir = sys.argv[2]
    out_dir = sys.argv[3]

    # cutoff = 75
    # indir = 'by_mag'
    # outdir = 'by_mag_filtered'

    
    filter_pathway(cutoff, mag_dir, out_dir)