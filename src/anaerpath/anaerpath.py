import requests
import os
import pandas
from Bio import SeqIO
import pandas as pd
import subprocess
import re
import glob

from .utils import *
from .scripts.KO_json_2_KO_table import transformJson2table, koTable2JsonIndex
from .scripts.compute_ko_expr import compute_ko_expr
from .scripts.filter_pathway import filter_pathway
from .scripts.search_pathway import search_pathway
from .scripts.pathway_sum import pathway_sum

class AnaerPath:
    '''
    Profile taxonomic genomes using a set of marker genes.
    '''
    def __init__(self, ko, gene, output, transcript=None, format='.fwd.fq/.rev.fq', process=30):
        self.ko = ko
        self.gene = gene
        self.transcript = transcript
        self.format = format
        self.process = process
        
        self.output = output
        self.path = os.path.dirname(os.path.realpath(__file__))

    def download_ko(self):
        '''
        Download ko00001 to data and make it plain.
        '''
        file = os.path.join(self.path, 'data', 'ko00001.json')
        try:
            with open(file, 'w') as w:
                for line in requests.get('https://www.kegg.jp/kegg-bin/download_htext?htext=ko00001.keg&format=json&filedir=').text:
                    w.write(line)

            if os.stat(file).st_size != 0:
                transformJson2table(file)
                koTable2JsonIndex(file)
                
        except:
            logger.warning('Fail to update ko00001, will use local version instead.')

    def merge_ko_with_pathway(self):
        '''
        Merge input ko file with updated pathway
        '''
        kegg_pathway_ko = pd.read_table(os.path.join(self.path, 'data', 'KEGG_pathway_ko.tsv'))
        ko = pd.merge(pd.read_table(self.ko), kegg_pathway_ko) # gene, ko
        
        

        ## split the ko table according to gene
        os.makedirs(os.path.join(self.output, '1.mag'), exist_ok=True)
        os.makedirs(os.path.join(self.output, '2.expr'), exist_ok=True)
        
        if self.transcript is not None:
            expression = pd.read_table(os.path.join(self.output, 'expression.tsv'))
            
        for file in glob.glob(os.path.join(self.gene, '*.genes')):
            filename = os.path.basename(file).rsplit('.genes', 1)[0]
            ids = set()
            with open(file) as handle:
                for record in SeqIO.parse(handle, 'fasta'):
                    ids.add(record.id)
            
            ko_subset = ko[ko.gene.isin(ids)].rename({'gene': 'gene_id', 'ko': 'ko'}, axis=1)
            ko_subset.to_csv(os.path.join(self.output, '1.mag', filename + '.txt'), sep='\t', index=False)

            ## create temp file for compute_ko_expr
            if self.transcript is None:
                ko_subset[['ko']].assign(count=1).to_csv(os.path.join(self.output, '2.expr', filename + '.txt'), sep='\t', index=False)
            else:
                ko_subset = pd.merge(ko_subset[['gene_id', 'ko']], expression).drop('gene_id', axis=1).to_csv(os.path.join(self.output, '2.expr', filename + '.txt'), sep='\t', index=False)
                
    def compute_ko_expr(self, completeness=75):
        os.makedirs(os.path.join(self.output, '3.pathway'), exist_ok=True)
        
        compute_ko_expr(
            os.path.join(self.output, '2.expr'),
            os.path.join(self.path, 'data', 'target_module.xlsx'),
            os.path.join(self.output, '3.pathway'))
                        
        filter_pathway(
            completeness,
            os.path.join(self.output, '3.pathway', 'by_mag'),
            os.path.join(self.output, '3.pathway', 'by_mag_filtered'))
            
        filter_pathway(
            completeness,
            os.path.join(self.output, '3.pathway', 'by_module'),
            os.path.join(self.output, '3.pathway', 'by_module_filtered'))
        
        os.makedirs(os.path.join(self.output, '4.pathsum'), exist_ok=True)
        search_pathway(
            os.path.join(self.output, '3.pathway', 'by_mag_filtered'),
            os.path.join(self.output, '2.expr'),
            os.path.join(self.output, '1.mag'),
            os.path.join(self.path, 'data', 'pathway_summary.xlsx'),
            os.path.join(self.output, '4.pathsum'))

        pathway_sum(
            os.path.join(self.output, '3.pathway', 'by_module_filtered'),
            self.output)


    def run_rsem(self, method):
        if self.transcript is None:
            pass
        else:
            ## cat all files
            os.makedirs(os.path.join(self.output, '0.rsem'), exist_ok=True)
            os.makedirs(os.path.join(self.output, '0.rsem'), exist_ok=True)
            with open(os.path.join(self.output, '0.rsem', 'genes.fa'), 'w') as outfile:
                for file in glob.glob(os.path.join(self.gene, '*.genes')):
                    with open(file) as infile:
                        for line in infile:
                            outfile.write(line)
        
            ## index
            subprocess.run([
                'rsem-prepare-reference',
                '--bowtie2',
                os.path.join(self.output, '0.rsem', 'genes.fa'),
                os.path.join(self.output, '0.rsem', 'genes.fa.rsem')
            ], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

            files = glob.glob(os.path.join(self.transcript, '*'))
            samples = {re.sub(self.format, '', x) for x in files}
            fwd, rev = self.format.split('|')

            for sample in samples:
                subprocess.run([
                    'rsem-calculate-expression',
                    '--paired-end',
                    '-p', str(self.process),
                    '--bowtie2',
                    '--no-bam-output',
                    sample + fwd,
                    sample + rev,
                    os.path.join(self.output, '0.rsem', 'genes.fa.rsem'),
                    os.path.join(self.output, '0.rsem', os.path.basename(sample))
                ], stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)

            expression = []
            for file in glob.glob(os.path.join(self.output, '0.rsem', '*.genes.results')):
                expression.append(pd.read_table(file)[['gene_id', method]].assign(sample = os.path.basename(file).split('.genes.results')[0]))

            expression = pd.concat(expression).set_index(['gene_id', 'sample']).unstack().fillna(0).droplevel(0, axis='columns')
            expression.to_csv(os.path.join(self.output, 'expression.tsv'), sep='\t')

    def run(self, completeness=75, skip_update=False, skip_rsem=False, method='TPM'):
        '''
        Run the pipeline.
        '''
        if not skip_rsem:
            logger.info('Running RSEM.')
            self.run_rsem(method=method)
        
        if not skip_update:
            logger.info('Updating ko00001.')
            self.download_ko()

        logger.info('Aggregating.')
        self.merge_ko_with_pathway()
        self.compute_ko_expr(completeness=completeness)
        
        logger.info('Done')