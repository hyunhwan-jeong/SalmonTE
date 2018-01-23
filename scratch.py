import os, sys, logging
from SalmonTE import *
from pathlib import Path
logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)


def get_ext(file_name):
    return os.path.basename(".".join(file_name.split('.')[1:]))


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
                name, anno = line[1:].strip().split("\t")[:-1]
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

build_salmon_index("reference/mm/mm.fa", "mm", True)
build_salmon_index("/Users/hwan/Downloads/dr.fa", "dr", True)