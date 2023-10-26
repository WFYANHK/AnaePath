import sys
import os

from argparse import ArgumentParser, SUPPRESS
from . import __version__
from .utils import logger
from .anaerpath import AnaerPath


def cli(argv=sys.argv):
    '''
    Entry point for command line interface.
    '''
    parser = ArgumentParser(description='anaerpath: ...', add_help=False)
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')

    required.add_argument(
        '-k',
        '--ko',
        metavar='FILE',
        required=True,
        help='File returned by GhostKOALA.')

    required.add_argument(
        '-g',
        '--gene',
        metavar='DIR',
        required=True,
        help='Folder containing genes (ORFs) predicted by prodigal.')

    required.add_argument(
        '-o',
        '--output',
        metavar='DIR',
        help='Output folder.')
    
    optional.add_argument(
        '-t',
        '--transcript',
        metavar='DIR',
        help='Folder containing metatranscriptomic files, must be paired, if given then use RSEM to compute expression. [None]')
    
    optional.add_argument(
        '-p',
        '--process',
        metavar='INT',
        type=int,
        default=30,
        help='Number of processes. [30]')
    
    optional.add_argument(
        '-f',
        '--format',
        metavar='STR',
        type=str,
        default='_fwd.fq|_rev.fq',
        help='Format of metatranscriptomic files, others can be _1.fq|_2.fq or _R1.fq|_R2.fq. [_fwd.fq|_rev.fq]')

    optional.add_argument(
        '-c',
        '--completeness',
        metavar='FLOAT',
        type=float,
        default=75,
        help='Pathway completeness cutoff. [75]')
    
    optional.add_argument(
        '-m',
        '--method',
        metavar='STR',
        type=str,
        default='TPM',
        help='Normalization method for RSEM output, can be RPKM or TPM. [TPM]')
    
    optional.add_argument(
        '--skip-update',
        action='store_true',
        help='Skip updating ko00001.')
    
    optional.add_argument(
        '--skip-rsem',
        action='store_true',
        help='Skip running RSEM.')

    parser.add_argument('-v', '--version', action='version', version=__version__, help=SUPPRESS)
    parser.add_argument('-h', '--help', action='help', help=SUPPRESS)

    opt = parser.parse_args(argv[1:])
    run(opt)


def run(opt):
    '''
    Sanity check of options.
    '''
    ## check for output folder
    if not os.path.isdir(opt.output):
        os.makedirs(opt.output, exist_ok=True)
    else:
        logger.warning('Folder <{}> exists. Files will be overwritten.'.format(opt.output))

    ## run
    AnaerPath(
        ko = opt.ko,
        gene = opt.gene,
        output = opt.output,
        transcript = opt.transcript,
        format = opt.format,
        process = opt.process
    ).run(completeness = opt.completeness, method=opt.method, skip_update = opt.skip_update, skip_rsem = opt.skip_rsem)

if __name__ == '__main__':
    cli(sys.argv)
