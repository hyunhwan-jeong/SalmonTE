#!/usr/bin/env python3
"""SalmonTE - Ultra-Fast and Scalable Quantification Pipeline of Transcript Abundances from Next Generation Sequencing Data

Usage:
    SalmonTE.py quant [--reference=genome] [--outpath=inpath] FILE...
    SalmonTE.py test [--pheno=type] [--mednorm] [--inpath=inpath] [--outpath=outpath]
    SalmonTE.py (-h | --help)
    SalmonTE.py --version

Options:
    -h --help     Show this screen.
    --version     Show version.
"""
from docopt import docopt


def SalmonTE(arguments):
    pass

if __name__ == '__main__':
    arguments = docopt(__doc__, version='SalmonTE 0.1')
    print(arguments)
    SalmonTE(arguments)
