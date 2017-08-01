#!/usr/bin/env python3
"""SalmonTE - Ultra-Fast and Scalable Quantification Pipeline of Transcript Abundances from Next Generation Sequencing Data

Usage:
    SalmonTE.py quant [--reference=genome] [--outpath=outpath] FILE...
    SalmonTE.py test [--pheno=type] [--mednorm] [--inpath=inpath] [--outpath=outpath]
    SalmonTE.py (-h | --help)
    SalmonTE.py --version

Options:
    -h --help     Show this screen.
    --version     Show version.
"""
from docopt import docopt
from glob import glob
from itertools import combinations
import logging
import sys
import os
import tempfile


def is_fastq(file_name):
    file_name = file_name.lower()
    if file_name.endswith("fastq"): return True
    if file_name.endswith("fq"): return True
    if file_name.endswith("fq.gz"): return True
    if file_name.endswith("fastq.gz"): return True
    return False


def get_basename_noext(file_name):
    return os.path.basename(file_name.split('.')[0])


def get_ext(file_name):
    return os.path.basename(".".join(file_name.split('.')[1:]))


def longest_prefix(a, b):
    a, b = (b,a) if len(a) > len(b) else (a,b)
    return b[:max([i for i in range(len(a),-1,-1) if b.startswith(a[:i])]+[0])]


def get_first_readid(file_name):
    with open(file_name, "r") as inp:
        return inp.readline().split()[0]


def collect_FASTQ_files(FILE):
    fastq_files = set()
    for file in FILE:
        collected_files = glob(file)
        if len(collected_files) == 0:
            logging.error("Failed to read {}".format(file))
            sys.exit(1)

        elif len(collected_files) == 1 and not is_fastq(collected_files[0]):
            collected_files = glob(collected_files[0]+"/*")

        for col_file in collected_files:
            if not is_fastq(col_file):
                logging.error("Failed to read {}".format(file))
                sys.exit(1)

        fastq_files |= set([os.path.abspath(file) for file in collected_files])

    tmp_dir = tempfile.mkdtemp()

    paired = dict([ (file,set()) for file in fastq_files ])
    for file_a, file_b in combinations(fastq_files, 2):
        if get_ext(file_a) != get_ext(file_b):
            logging.info("File extensions of all files must be same.")
            sys.exit(1)

        if get_first_readid(file_a) == get_first_readid(file_b):
            paired[file_a].add(file_b)
            paired[file_b].add(file_a)

    is_paired = True
    for file in paired:
        if len(paired[file]) == 1: is_paired = True
        elif len(paired[file]) == 0: is_paired = False
        elif is_paired and len(paired[file]) > 1:
            logging.error("One of files can be paired to multiple files")
            sys.exit(1)
        elif not is_paired and len(paired[file]) > 0:
            logging.error("The pair should not exists")
            sys.exit(1)

    if is_paired:
        for file, matched in paired.items():
            if len(matched) == 0:
                logging.error("Input dataset is supposed as a paired-ends dataset, but some files do not have thier pair.")
                sys.exit(1)

        logging.info("The input dataset is considered as a paired-ends dataset.")
        for file in sorted(paired):
            a = file
            b = list(paired[file])[0]
            if a > b: continue
            trim_a = "_".join(get_basename_noext(a).split("_")[:-1]) + "_R1." + ".".join(os.path.basename(a).split('.')[1:])
            trim_b = "_".join(get_basename_noext(b).split("_")[:-1]) + "_R2." + ".".join(os.path.basename(b).split('.')[1:])
            os.symlink(os.path.abspath(a), os.path.join(tmp_dir, trim_a))
            os.symlink(os.path.abspath(b), os.path.join(tmp_dir, trim_b))
    else:
        for file in sorted(fastq_files):
            logging.info("The input dataset is considered as a single-end dataset.")
            os.symlink(os.path.abspath(file), os.path.join(tmp_dir, get_basename_noext(file)))

    ret = dict()
    ret["paired"] = is_paired
    ret["inpath"] = tmp_dir
    ret["num_fastq"] = int(len(fastq_files) / (is_paired*2)) if is_paired else len(fastq_files)
    return ret

def run(args):
    if args['quant']:
        print(args)
        logging.info("Starting quantification mode")
        logging.info("Collecting FASTQ files...")
        inp_param = collect_FASTQ_files(args['FILE'])
        logging.info("Collectd {} FASTQ files.".format(inp_param["num_fastq"]))
        logging.info("Quantification has been finished.")

if __name__ == '__main__':
    logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)

    args = docopt(__doc__, version='SalmonTE 0.1')
    run(args)
