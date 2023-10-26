#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import glob
import os

import numpy as np
import re
import sys

def _parse(string):
    '''
    Parse substring, return a float
    '''
    ## if submodule inside parenthesis, replace it with chain
    if ' ' in string:
        string = string.replace(' ', '+')
        
    ## for each parts, mean, ignore missing parts
    parts = re.sub(r'\(|\)', '', string).split(',')
    
    arrays = []
    for part in parts:
        array = np.array(re.split(r'\+|\-', part), dtype=float)
        
        if np.isnan(array).all():
            arrays.append(np.nan)
        else:
            arrays.append(np.nanmean(array))
    
    ## for all parts, sum, ignore missing parts
    if np.isnan(arrays).all():
        return(np.nan)
    else:
        return(np.nansum(arrays))
    
def compute_ko_expr(mag_dir, module_file, out_dir):
    '''
    Compute expression for each module/mag/sample combination
    '''
    mags = glob.glob(os.path.join(mag_dir, '*.txt'))
    modules = pd.read_excel(module_file)

    results = []
    for mag in mags:
        ko_sample_expr = pd.read_csv(mag, sep='\t').groupby(['ko']).sum()
        # print(ko_sample_expr)
        for sample in ko_sample_expr.columns:
            for module, definition in zip(modules['Module'], modules['DEFINITION']):
                ko_expr = ko_sample_expr[sample].to_dict()

                if '--' in definition or '\n' in definition:
                    print('Module <{}> contains merging or branching events, skip'.format(module))
                    continue

                ## replace ko with expr, fill missing with nan
                for ko in re.findall(r'K\d{5}', definition):
                    definition = definition.replace(ko, str(ko_expr.get(ko, np.nan)))
                # print(definition)

                ## do calculation for each parenthesis
                while '(' in definition:
                    for part in re.findall('\([^()]+\)', definition):
                        definition = definition.replace(part,  str(_parse(part)))

                ## for all submodules, mean, but don't ignore missing submodules
                submodules = definition.split(' ')
                # print(mag, module, sample, definition)
                array = np.array([_parse(submodule) for submodule in submodules], dtype=float)
                results.append([mag, module, sample, np.isnan(array).sum(), len(array), np.nansum(array)/len(array)])


    df = pd.DataFrame(results, columns = ['MAG', 'Pathway_name', 'sample', 'Steps_lacking', 'Steps_needed', 'expression'])
    df['MAG'] = df['MAG'].str.replace('\.txt$|.*/', '', regex=True)
    df['Pathway_completeness'] = round(100 * (1 - df['Steps_lacking']/df['Steps_needed']) , 2)
    df['expression'] = round(df['expression'], 2)

    ## save files
    for module, temp in pd.DataFrame(df).groupby('Pathway_name'):
        os.makedirs(os.path.join(out_dir, 'by_module'), exist_ok=True)
        out = os.path.join(out_dir, 'by_module', module + '.xlsx')
        sheet = temp.pivot(index=['MAG','Steps_lacking','Steps_needed','Pathway_completeness'], columns='sample', values='expression')
        # sheet.reindex(columns = sorted(sheet.columns, key=lambda x: -float(x[1:].split('_')[0]) + 0.1 * float(x.split('_')[-1]))).to_excel(out)
        sheet.to_excel(out)

    for mag, temp in pd.DataFrame(df).groupby('MAG'):
        os.makedirs(os.path.join(out_dir, 'by_mag'), exist_ok=True)
        out = os.path.join(out_dir, 'by_mag', mag + '.xlsx')
        sheet = temp.pivot(index=['Pathway_name','Steps_lacking','Steps_needed','Pathway_completeness'], columns='sample', values='expression')
        sheet.to_excel(out)
        # sheet.reindex(columns = sorted(sheet.columns, key=lambda x: -float(x[1:].split('_')[0]) + 0.1 * float(x.split('_')[-1]))).to_excel(out)

    out = os.path.join(out_dir, 'all_MAG.xlsx')
    sheet = pd.DataFrame(df).pivot(index=['Pathway_name','MAG','Steps_lacking','Steps_needed','Pathway_completeness'], 
                           columns='sample', values='expression')
    # sheet.reindex(columns = sorted(sheet.columns, key=lambda x: -float(x[1:].split('_')[0]) + 0.1 * float(x.split('_')[-1]))).to_excel(out, merge_cells=False)
    sheet.to_excel(out, merge_cells=False)

if __name__ == '__main__':
    mag_dir = sys.argv[1]
    module_file = sys.argv[2]
    out_dir = sys.argv[3]

    # mag_dir = '1.MAG_file'
    # module_file = '2.target_module.xlsx'
    # out_dir = 'test'

    os.makedirs(out_dir, exist_ok=True)
    compute_ko_expr(mag_dir, module_file, out_dir)
