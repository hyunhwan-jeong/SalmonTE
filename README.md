## Change Logs
* Version 0.2: Support two different type of expressions with `--exprtype` option 

## What is SalmonTE?
`SalmonTE` is an ultra-Fast and Scalable Quantification Pipeline of Transpose Element (TE) Abundances from Next Generation Sequencing Data. It comes with [Salmon](https://github.com/COMBINE-lab/salmon) which is a fast and accurate transcriptome quantification method. You can read the details of the pipeline and an example of real data study in [my recent publisehd paper in PSB 2018](http://www.worldscientific.com/doi/10.1142/9789813235533_0016).

## Which data you need to run `SalmonTE`? Why I have to use it?
* You only need to have a set of FASTQ files and phenotype data. Furthermore, **`SalmonTE` automatically decided wether your dataset is paired-ends reads or not.** 
* phenotype can be a numeric data or a categorical data. Based on the data type of the phenotype, **`SalmonTE` will run differential expression analysis if user's input contains 'control' or 'case' string as the phenotype**. Otherwise, `SalmonTE` will run `LM` for the data.
* Unlikely other TE analysis tools, `SalmonTE` gives you various visualized output. It must be helpful to your research.

## Requirements & Installation
To use `SalmonTE` `python` and `R` must be installed before running it.

* **Note**: Currently, running `SalmonTE` on MacOS has an issue, and we are try to fix it soon. Thus, we recommend to use `linux` environment to play it.

* Install `python` and `R` packages

For `python`:

Run following line in your console
```
pip3 install snakemake docopt --user
```

For `R`:
Run following lines in `R` console.

```
install.packages(c("tidyverse", "cowplot", "scales", "WriteXLS"))
source("https://bioconductor.org/biocLite.R")
biocLite(c("DESeq2", "tximport"))
```


* Clone the repository 

```
git clone https://github.com/hyunhwaj/SalmonTE
```

* Add `PATH` of SalmonTE to your `.bashrc` file:

```
export PATH=$PATH:/PATH_OF_SALMON_TE/
```

* Re log-in to terminal or use `source` command:

```
source ~/.bashrc
```

## How to use it?

```
Usage:
    SalmonTE.py quant [--reference=genome] [--outpath=outpath] [--num_threads=numthreads] [--exprtype=exprtype] FILE...
    SalmonTE.py test [--inpath=inpath] [--outpath=outpath] [--tabletype=tabletype] [--figtype=figtype]
    SalmonTE.py (-h | --help)
    SalmonTE.py --version

Options:
    -h --help     Show this screen.
    --version     Show version.
```

## An real example of SalmonTE usage with command line 
```
SalmonTE.py quant --reference=hs example
# Prior to run statistical analysis, open `phenotype.csv` and put your phenotype/covariate data.
vi SalmonTE_output/phenotype.csv 
SalmonTE.py --inpath=SalmonTE_output --outpath=tmp --tabletype=csv --figtype=png
```
