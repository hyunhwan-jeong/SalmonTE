# What is SalmonTE?

SalmonTE is an ultra-Fast and Scalable Quantification Pipeline of Transpose Element (TE) Abundances from Next Generation Sequencing Data

# How to use it?

```
Usage:
    SalmonTE.py quant [--reference=genome] [--outpath=outpath] [--num_threads=numthreads] FILE...
    SalmonTE.py test [--inpath=inpath] [--outpath=outpath] [--tabletype=tabletype] [--figtype=figtype]
    SalmonTE.py (-h | --help)
    SalmonTE.py --version

Options:
    -h --help     Show this screen.
    --version     Show version.
```

# What you need to run SalmonTE?

* You only need to have a set of FASTQ files and phenotype data
* phenotype can be a numeric data or a categorical data

# An real example of SalmonTE usage with command line 

```
SalmonTE.py quant quant --reference=hs example
SalmonTE.py --inpath=SalmonTE_output --outpath=tmp --tabletype=csv --figtype=png
```
