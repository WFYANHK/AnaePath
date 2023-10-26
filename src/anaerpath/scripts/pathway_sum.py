#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import glob
import os

def pathway_sum(indir, outdir):
    if not os.path.isdir(outdir):
        os.makedirs(outdir, exist_ok=True)

    tmp = []
    for file in glob.glob(os.path.join(indir, '*.xlsx')):
        filename = os.path.basename(file)
        # print('Processing ' + filename + '...')
        df = pd.read_excel(file)
        df['Pathway_name'] = filename
        tmp.append(df)

    a = pd.concat(tmp).drop(['Steps_lacking', 'Steps_needed', 'Pathway_completeness','MAG'], axis=1).groupby('Pathway_name').sum()
    b = pd.concat(tmp).groupby('Pathway_name')['MAG'].apply(lambda x: ','.join(x))
    pd.merge(b,a, left_index=True, right_index=True).to_csv(os.path.join(outdir, 'pathway.tsv'), sep='\t')