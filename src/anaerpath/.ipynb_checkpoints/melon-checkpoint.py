import glob
import pandas as pd
from collections import defaultdict

from .utils import *

class GenomeProfiler:
    '''
    Profile taxonomic genomes using a customized set of marker-genes.
    '''
    def __init__(self, file, output='tmp', threads=32):
        self.file = file
        self.output = output
        self.threads = threads

        self.bset = {'l11', 'l2', 'l13', 'l19', 'l20', 's7', 's15', 's6'}
        self.aset = {'l11', 'l2', 'l15e', 'l16/l10e', 'l44e', 's19e', 's3ae', 's9'}
        self.nset = set()

        self.ranks = ['Superkingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']

        ## genome copies
        self.copies = {'bacteria': 0, 'archaea': 0}

        ## temporary variables for diamond hits or minimap mappings
        self.hits, self.maps = [], []

        ## taxonomy assignments
        self.assignments = {}


    def run_kraken(self, db_kraken):
        '''
        Run kraken2 for pre-filtering.
        '''
        subprocess.run([
            'kraken2', 
            '--threads', str(self.threads), 
            '--report', get_filename(self.file, self.output, '.rep.tmp'), 
            '--output', get_filename(self.file, self.output, '.out.tmp'),
            '--db', db_kraken, 
            self.file,
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


    def parse_kraken(self):
        '''
        Parse kraken2's output by 
            1: recording the ids of eukaryota (2759), viruses (10239), other entries (2787854).
            2: adding the ids of the sequences to a negative search set. 
        '''
        seqid, taxid = set(), set()
        record = False
        with open(get_filename(self.file, self.output, '.rep.tmp')) as f:
            for line in f:
                ls = line.rstrip().split('\t')
                if ls[3] in {'D', 'R1'}:
                    record = ls[5].strip() in {'Eukaryota', 'Viruses', 'other entries'}

                if record:
                    taxid.add(ls[4])

        with open(get_filename(self.file, self.output, '.out.tmp')) as f:
            for line in f:
                ls = line.rstrip().split('\t')
                if ls[2] in taxid:
                    seqid.add(ls[1])

        self.nset.update(seqid)


    def run_diamond(self, db):
        '''
        Run diamond to get total genome copies.
        '''
        outfmt = ['qseqid', 'sseqid', 'pident', 'length', 'qlen', 'qstart', 'qend', 'slen', 'sstart', 'send', 'evalue', 'bitscore']
        subprocess.run([
            'diamond', 'blastx',
            '--db', os.path.join(db, 'prot.fa'),
            '--query', self.file,
            '--out', get_filename(self.file, self.output, '.hit.tmp'),
            '--outfmt', '6', *outfmt,
            '--evalue', '1e-10',
            '--subject-cover', '75',
            '--threads', str(self.threads),
            '--range-culling', '-F', '15',
            '--max-hsps', '0',
            '--range-cover', '25'
        ], check=True)


    def parse_diamond(self):
        '''
        Parse diamond's output and record the hits.
        '''
        qrange = defaultdict(set)
        srange = {}

        with open(get_filename(self.file, self.output, '.hit.tmp')) as f:
            for line in f:
                ls = line.rstrip().split('\t')
                qseqid, sseqid = ls[0], ls[1]
                qstart, qend = sort_coordinate(int(ls[5]), int(ls[6]))
                sstart, send = sort_coordinate(int(ls[8]), int(ls[9]))
                slen = int(ls[7])

                ## bypass non-prokaryotic reads
                if qseqid not in self.nset:
                    if (
                        qseqid not in qrange or
                        all(compute_overlap((qstart, qend, *x), max) < 0.25 for x in qrange.get(qseqid))
                    ):
                        qrange[qseqid].add((qstart, qend))

                        ss = sseqid.split('-')
                        gene, kingdom = ss[0], ss[-2]
                        if (
                            (gene in self.bset and kingdom == 'bacteria') or
                            (gene in self.aset and kingdom == 'archaea')
                        ):
                            if sseqid not in srange:
                                srange[sseqid] = np.zeros(slen)
                            srange[sseqid][range(sstart, send)] += 1

                            ## append qseqid and coordinates for back-tracing
                            self.hits.append([qseqid, kingdom, gene, qstart, qend])

        for key, val in srange.items():
            cut = round(len(val)/4) # in this case same as counting hits since all hits have subject cover > 75%
            kingdom = key.split('-')[-2]
            if kingdom == 'bacteria':
                self.copies[kingdom] += np.mean(val[cut:-cut])/len(self.bset)
            else:
                self.copies[kingdom] += np.mean(val[cut:-cut])/len(self.aset)


    def run_minimap(self, db):
        '''
        Run minimap2 to get taxonomic profiles.
        '''
        with open(get_filename(self.file, self.output, '.seq.tmp'), 'w') as w:
            w.write(extract_sequences(self.file, {x[0] for x in self.hits}))

        genes = defaultdict(set)
        for hit in self.hits:
            genes[hit[1] + '.' + hit[2].replace('/', '_')].add(hit[0])

        with open(get_filename(self.file, self.output, '.map.tmp'), 'w') as w:
            for key, val in genes.items():
                sequences = extract_sequences(get_filename(self.file, self.output, '.seq.tmp'), val)
                subprocess.run([
                    'minimap2',
                    '-cx', 'map-ont',
                    '-t', str(self.threads),
                    '-N', '1000',
                    '-f', '0',
                    os.path.join(db, 'nucl.' + key + '.fa'), '-',
                ], check=True, stdout=w, stderr=subprocess.DEVNULL, input=sequences, text=True)


    def parse_minimap(self):
        '''
        Parse minimap2's output and record the mappings.
        '''
        coordinates = defaultdict(set)
        for i in self.hits:
            coordinates[i[0]].add(tuple(i[-2:]))

        with open(get_filename(self.file, self.output, '.map.tmp')) as f:
            for line in f:
                ls = line.split('\t')
                qstart, qend, qseqid, sseqid = int(ls[2]), int(ls[3]), ls[0], ls[5]
                accession = sseqid.rsplit('_',3)[0]

                AS = float(ls[14].split(':')[-1]) # alignment score
                for i in coordinates.get(qseqid):
                    self.maps.append([qseqid, sseqid, accession, AS, qstart, qend] + list(i))

        ## filter out non-overlapping mappings
        self.maps = [x for x in self.maps if compute_overlap(x[-4:])>0]


    def postprocess(self, db):
        '''
        Post-processing and label reassignment using EM.
        '''
        lineage = {}
        with open(os.path.join(db, 'metadata.tsv')) as f:
            next(f)
            for line in f:
                ls = line.rstrip().split('\t')
                lineage[ls[0]] = ';'.join(ls[1:])

        columns = ['qseqid', 'sseqid', 'accession', 'AS', 'qstart', 'qend', 'pstart', 'pend']
        data = pd.DataFrame(self.maps, columns=columns).sort_values(['qseqid', 'AS'], ascending=[False, False])
        data['lineage'] = data['accession'].map(lineage)

        assert (
            data.lineage.isnull().sum()==0
        ), 'Not all accessions have a valid taxonomy.'

        ## create a score matrix
        data = data.groupby(['qseqid', 'lineage'], as_index=False).first()
        data['ratio'] = data['AS'] / data.groupby('qseqid')['AS'].transform('max')
        data['score'] = pd.cut(data['ratio'], bins=[0.9, 0.99, 0.999, 1], labels=[1/4, 1/2, 1], right=True).astype(float)
        matrix = data[data['score'].notnull()].set_index(['qseqid', 'lineage'])['score'].unstack().fillna(0)

        ## run EM
        assignments = reassign_taxonomy(matrix.to_numpy())

        ties = defaultdict(list)
        for i, j in zip(matrix.index, assignments):
            assignment = [matrix.columns[x] for x in j]
            if len(assignment)>1:
                ties[tuple(assignment)].append(i)
            else:
                self.assignments[i] = assignment[0]

        ## resolve equal probability cases using alignment score
        for key, val in ties.items():
            target = data[(data['qseqid'].isin(val)) & (data['lineage'].isin(key))].groupby('lineage')['AS'].mean()
            target = target.sort_values(ascending=False).index[0] # select a single one, though may still tie
            for i in val:
                self.assignments[i] = target

    def run(self, db, db_kraken, skip_profile=False, skip_clean=False):
        '''
        Run the pipeline.
        '''
        if db_kraken is not None:
            logger.info('Filtering reads ...')
            self.run_kraken(db_kraken)
            self.parse_kraken()
            logger.info('... removed {} putative non-prokaryotic reads.'.format(len(self.nset)))

        logger.info('Estimating genome copies ...')
        self.run_diamond(db)
        self.parse_diamond()
        logger.info('... found {} genome copies (bacteria: {}; archaea: {}).'.format(
            sum(self.copies.values()), self.copies.get('bacteria', 0), self.copies.get('archaea', 0)))

        if not skip_profile:
            logger.info('Assigning taxonomy ...')
            self.run_minimap(db)
            self.parse_minimap()
            self.postprocess(db)

            ## generate a profile output
            columns = ['qseqid', 'kingdom', 'lineage']
            profile = pd.DataFrame([(x[0], x[1], self.assignments.get(x[0])) for x in self.hits], columns=columns)

            ## fill missing according to hits
            for kingdom in ['bacteria', 'archaea']:
                if kingdom == 'bacteria':
                    replacement = ';'.join(['2|Bacteria'] + ['0|unclassified Bacteria ' + x.lower() for x in self.ranks[1:]])
                else:
                    replacement = ';'.join(['2157|Archaea'] + ['0|unclassified Archaea ' + x.lower() for x in self.ranks[1:]])
                profile[(profile['lineage'].isnull()) & (profile['kingdom'] == kingdom)] = replacement

            ## split lineage to seven plain taxonomic ranks
            profile[self.ranks] = profile['lineage'].str.split(';', expand=True)

            ## aggregate to species-level genome copies and relative abundance
            profile = profile.groupby(self.ranks, as_index=False).size()
            profile['Copies'] = sum(self.copies.values()) * profile['size'] / profile['size'].sum()
            profile['Abundance'] = profile['Copies'] / profile['Copies'].sum()

            profile = profile.sort_values(['Copies'] + self.ranks)[self.ranks + ['Copies', 'Abundance']]
            profile.to_csv(get_filename(self.file, self.output, '.tsv'), sep='\t', index=False)

            ## check number of species obtained
            species = profile[~profile['Species'].str.contains('0|unclassified', regex=False)]
            species = species.groupby(['Superkingdom']).size().to_dict()
            logger.info('... identified {} unique species (bacteria: {}; archaea: {}).'.format(
                sum(species.values()), species.get('2|Bacteria', 0), species.get('2157|Archaea', 0)))

        else:
            profile = pd.DataFrame([
                ['2|Bacteria', self.copies.get('bacteria')],
                ['2157|Archaea', self.copies.get('archaea')]], columns = ['Superkingdom', 'Copies'])

            profile['Abundance'] = profile['Copies'] / profile['Copies'].sum()
            profile.sort_values('Copies').to_csv(get_filename(self.file, self.output, '.tsv'), sep='\t', index=False)

        ## clean up
        if not skip_clean:
            for f in glob.glob(get_filename(self.file, self.output, '.*.tmp')):
                os.remove(f)
