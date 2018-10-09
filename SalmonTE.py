#!/usr/bin/env python3
"""SalmonTE - Ultra-Fast and Scalable Quantification Pipeline of Transcript Abundances from Next Generation Sequencing Data

Usage:
    SalmonTE.py index [--ref_name=ref_name] (--input_fasta=fa_file) [--te_only]
    SalmonTE.py quant [--reference=genome] [--outpath=outpath] [--num_threads=numthreads] [--exprtype=exprtype] FILE...
    SalmonTE.py test [--inpath=inpath] [--outpath=outpath] [--tabletype=tabletype] [--figtype=figtype]
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
import gzip

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
    return os.path.basename(".".join(file_name.split('.')[-1:]))


def longest_prefix(a, b):
    a, b = (b,a) if len(a) > len(b) else (a,b)
    return b[:max([i for i in range(len(a),-1,-1) if b.startswith(a[:i])]+[0])]


def get_first_readid(file_name):
    if file_name.lower().endswith(".gz"):
        with gzip.open(file_name, "rb") as inp:
            return inp.readline().decode("utf8").replace("/"," ").split()[0]
    else:
        with open(file_name, "r") as inp:
            return inp.readline().replace("/"," ").split()[0]


def correct_ext(ext):
    return "fastq.gz" if "gz" in ext else "fastq"

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
        """
        if get_ext(file_a) != get_ext(file_b):
            logging.info("File extensions of all files must be same.")
            sys.exit(1)
        """

        if get_first_readid(file_a) == get_first_readid(file_b):
            paired[file_a].add(file_b)
            paired[file_b].add(file_a)

    is_paired = True
    for file in paired:
        if len(paired[file]) == 0: is_paired = False



    for file in paired:
        if is_paired and len(paired[file]) > 1:
            logging.error("One of input files can be paired to multiple files")
            sys.exit(1)
        elif not is_paired and len(paired[file]) > 0:
            logging.error("A paired-end sample and a single-end sample are placed together.")
            sys.exit(1)

    file_list = []
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
            trim_a = "_".join(get_basename_noext(a).split("_")[:-1]) + "_R1.{}".format(correct_ext(os.path.basename(a).split('.')[1:]))
            trim_b = "_".join(get_basename_noext(a).split("_")[:-1]) + "_R2.{}".format(correct_ext(os.path.basename(b).split('.')[1:]))
            os.symlink(os.path.abspath(a), os.path.join(tmp_dir, trim_a))
            os.symlink(os.path.abspath(b), os.path.join(tmp_dir, trim_b))
            file_list.append([(trim_a, trim_b)])
    else:
        for file in sorted(fastq_files):
            logging.info("The input dataset is considered as a single-end dataset.")
            file_name = os.path.join(tmp_dir, get_basename_noext(file)) + ".{}".format(correct_ext(os.path.basename(file).split('.')[1:]))
            os.symlink(os.path.abspath(file), file_name)
            file_list.append(os.path.basename(file))

    ret = dict()
    ret["paired"] = is_paired
    ret["inpath"] = tmp_dir
    ret["num_fastq"] = int(len(fastq_files) / (is_paired*2)) if is_paired else len(fastq_files)
    ret["file_list"] = file_list
    return ret

def run_salmon(param):
    import snakemake
    snakefile = os.path.join(os.path.dirname(__file__), "snakemake/Snakefile.paired" if param["paired"] else "snakemake/Snakefile.single")

    snakemake.snakemake(
        snakefile=snakefile,
        config={
            "input_path": param["inpath"],
            "output_path": param["--outpath"],
            "index": param["--reference"],
            "salmon": os.path.join(os.path.dirname(__file__),"salmon/{}/bin/salmon"),
            "num_threads" : param["--num_threads"],
            "exprtype": param["--exprtype"]
        }
    )

    with open(os.path.join(param["--outpath"], "EXPR.csv" ), "r") as inp:
        sample_ids = inp.readline().strip().split(',')[1:]
    with open(os.path.join(param["--outpath"], "phenotype.csv" ), "w") as oup:
        oup.write("SampleID,phenotype\n")
        oup.write(
            "\n".join([s+","+"NA" for s in sample_ids]) + "\n"
        )


def build_salmon_index(in_fasta, ref_id, te_only):
    logging.info("Building Salmon Index")
    if not os.path.exists(in_fasta):
        logging.error("Input file does not exist.")
        sys.exit(1)
    ext = get_ext(in_fasta)
    if not ext.lower() in [ "fa", "fasta" ]:
        logging.error("Input file is not a FASTA file.")
        sys.exit(1)
    ref_path = os.path.join(os.path.dirname(__file__), "reference")
    out_path = os.path.join(ref_path, ref_id)

    if not os.path.exists(out_path):
        os.mkdir(out_path)

    clade_dict = {}
    with open(os.path.join(ref_path, "clades_extended.csv"), "r") as inp:
        for line in inp:
            anno, clade, repeat_class, category = line.strip().split(",")
            if te_only and category != "Transposable Element": continue
            clade_dict[anno] = [clade, repeat_class]


    file_fa = os.path.join(out_path, "ref.fa")
    file_clade = os.path.join(out_path, "clades.csv")

    with open(in_fasta, "r") as inp, \
         open(file_fa, "w") as oup_fa, \
         open(file_clade, "w") as oup_clade:
        oup_clade.write("name,class,clade\n")
        for line in inp:
            if line.startswith(">"):
                name, anno = line[1:].strip().split("\t")[:2]
                if anno in clade_dict:
                    clade, repeat_class = clade_dict[anno]
                    oup_fa.write(">{}\n".format(name))
                    oup_clade.write("{},{},{}\n".format(name,
                                                       clade,
                                                       repeat_class))
            else:
                if anno in clade_dict:
                    oup_fa.write(line.strip()+"\n")

    salmon_bin = os.path.join(os.path.dirname(__file__), "salmon/{}/bin/salmon").format(sys.platform)
    cmd = "{} index -t {} -i {} --type=quasi".format(salmon_bin, file_fa, out_path)
    os.system(cmd)
    logging.info("Building '{}' index was finished!".format(ref_id))


def run(args):
    if args['quant']:
        if args['--exprtype'] is None:
            args['--exprtype'] = "TPM"
        if args['--num_threads'] is None:
            args['--num_threads'] = 4
        if args['--outpath'] is None:
            args['--outpath'] = os.path.join(os.getcwd(), "SalmonTE_output/")
        if args['--reference'] is None or args['--reference'] == "hs":
            args['--reference'] = os.path.join(os.path.dirname(__file__), "reference/hs")
        elif os.path.exists(os.path.join(os.path.dirname(__file__), "reference", args['--reference'],"versionInfo.json")):
            args['--reference'] = os.path.join(os.path.dirname(__file__), "reference", args['--reference'])
        else:
            logging.error("Reference file is not found!")
            sys.exit(1)


        logging.info("Starting quantification mode")
        logging.info("Collecting FASTQ files...")
        param = {**args, **collect_FASTQ_files(args['FILE'])}
        logging.info("Collected {} FASTQ files.".format(param["num_fastq"]))
        logging.info("Quantification has been finished.")
        logging.info("Running Salmon using Snakemake")

        run_salmon(param)
        os.system("cp {}/clades.csv {}".format(args['--reference'],
                                               args['--outpath']))

    if args['index']:
        build_salmon_index(args['--input_fasta'], args['--ref_name'], args['--te_only'])

    if args['test']:
        import snakemake
        if args['--inpath'] is None:
            logging.error("Input path must be specified!")
            sys.exit(1)
        elif not os.path.exists(os.path.join(args['--inpath'], "EXPR.csv")):
            logging.error("Input path is specified incorrectly!")
            sys.exit(1)

        if args['--outpath'] is None:
            snakemake.utils.makedirs(os.path.join(os.getcwd(), "SalmonTE_output"))
        elif not os.path.exists(args['--outpath']):
            snakemake.utils.makedirs(args['--outpath'])

        if args['--tabletype'] is None:
            args['--tabletype'] = "xls"

        if args['--figtype'] is None:
            args['--figtype'] = "pdf"


        #  Rscript SalmonTE_Stats.R SalmonTE_output xls PDF tmp
        os.system("Rscript {} {} {} {} {}".format( os.path.join(os.path.dirname(__file__), "SalmonTE_Stats.R"),
                                        args["--inpath"],
                                        args['--tabletype'],
                                        args['--figtype'],
                                        args['--outpath']))




if __name__ == '__main__':
    logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)
    args = docopt(__doc__, version='SalmonTE 0.3')
    run(args)
