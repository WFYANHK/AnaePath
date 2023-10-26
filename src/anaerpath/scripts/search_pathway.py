#!/usr/bin/env python
# coding: utf-8

import glob
import sys
import os
import pandas as pd

def search_pathway(mag_dir, tpm_dir, ko_dir, path_file, out_dir):
    '''
    Combine filtered mag and ko - tpm
    '''
    ps = pd.read_excel(path_file)
    ps['Pathway_name'] = ps['Pathway_name'].str.rstrip()

    for file in glob.glob(os.path.join(mag_dir, '*.xlsx')):

        if '$' not in file: # skip temp. excel files
            mag = '.'.join(os.path.basename(file).split('.')[:-1])
            ## step 1 -> merge filtered mag with path_file, get ko_list
            tmp = pd.read_excel(file)
            tmp['Pathway_name'] = tmp['Pathway_name'].str.rstrip()
            ko_list = pd.merge(tmp, ps, on='Pathway_name')[['Pathway_name', 'KO', 'Step','Gene name', 'Enzyme name', 'EC number', 'Substrate', 'Product']]

            try:
                ## step 2 -> merge ko_list with ko, get gene id
                ko = pd.read_csv(os.path.join(ko_dir, mag+'.txt'), sep='\t')[['gene_id', 'ko']].rename({'ko':'KO', 'gene_id':'GeneID'}, axis=1)
                ko_list = pd.merge(ko_list, ko, on='KO')

                ## step 3 -> use gene id to get tpm
                tpm = pd.read_csv(os.path.join(tpm_dir, mag+'.txt'), sep='\t').rename({'ko':'KO'}, axis=1)
                ko_list = pd.merge(ko_list, tpm, on='KO')

                ## step 4 -> reorder columns & save
                ko_list = ko_list[['Pathway_name', 'Step', 'GeneID', 'KO', 'Gene name', 'Enzyme name', 'EC number', 'Substrate', 'Product'] + list(ko_list.columns)[9:]]
                ko_list.sort_values(['Pathway_name', 'Step']).to_excel(os.path.join(out_dir, mag + '.xlsx'), index=False)

            except Exception as e:
                print(e)

if __name__ == '__main__':

    mag_dir = sys.argv[1]
    tpm_dir = sys.argv[2]
    ko_dir = sys.argv[3]
    path_file = sys.argv[4]
    out_dir = sys.argv[5]

    # mag_dir = 'by_mag_filtered'
    # tpm_dir = 'input1'
    # ko_dir = 'input2'
    # path_file = '4.pathway_summary.xlsx'
    # out_dir = 'test'
    
    os.makedirs(out_dir, exist_ok=True)
    search_pathway(mag_dir, tpm_dir, ko_dir, path_file, out_dir)
