## Change Logs

* January 23, 2017: Added references of *Mus musculus*(mm) and *Danio rerio*(dr), added a function users to allow to build a customized index (`index`) mode, and fixed a minor bug.

* December 25, 2017: Fixed a bug for single-end reads dataset.

* December 14, 2017: Fixed issue of the `macOS`. We have figured out there is a problem of old version of `snakemake`. If you already install the package, and the version is not above `4.0.0` (You can check it with `snakemake --version`) then please update version with below command:
`pip3 install snakemake --user --upgrade`

* November 27, 2017: source code of PSB manuscript is out now - our [manuscript](https://github.com/hyunhwaj/SalmonTE-manuscript/) and [response letter](https://github.com/hyunhwaj/SalmonTE-response).

* November 21, 2017: Version 0.2 is out, `SalmonTE` now supports two different type of expressions with `--exprtype` option. Use `--exprtype=TPM` if you want to use TPM values for the statistical analysis. If you want to run differential expression analysis, then I highly recommend to use `--exprtype=count`. [Here is the nice answer why](https://support.bioconductor.org/p/98820/).

## What is SalmonTE?
`SalmonTE` is an ultra-Fast and Scalable Quantification Pipeline of Transpose Element (TE) Abundances from Next Generation Sequencing Data. It comes with [Salmon](https://github.com/COMBINE-lab/salmon) which is a fast and accurate transcriptome quantification method. You can read the details of the pipeline and an example of real data study in [my recent publisehd paper in PSB 2018](http://www.worldscientific.com/doi/10.1142/9789813235533_0016).

## Which data you need to run `SalmonTE`? Why I have to use it?
* You only need to have a set of FASTQ files and phenotype data. Furthermore, **`SalmonTE` automatically decided wether your dataset is paired-ends reads or not.** 
* phenotype can be a numeric data or a categorical data. Based on the data type of the phenotype, **`SalmonTE` will run differential expression analysis if user's input contains 'control' or 'case' string as the phenotype**. Otherwise, `SalmonTE` will run `LM` for the data.
* Unlikely other TE analysis tools, `SalmonTE` gives you various visualized output. It must be helpful to your research.

## Requirements & Installation
To use `SalmonTE` `python` and `R` must be installed before running it.

~* **Note**: Currently, running `SalmonTE` on MacOS has an issue, and we are try to fix it soon. Thus, we recommend to use `linux` environment to play it.~

* Install `python` and `R` packages

For `python`:

Run following line in your console
```
pip3 install snakemake docopt pandas --user
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

### Troubleshooting

**Q. I am using `SalmonTE` on `macOS` and `salmonTE` fails to run on `quant` mode with error messages:**

```
CalledProcessError in line xx of SOME_PATH:
Command ' set -euo pipefail;  ROOT_OF_SALMON_TE/SalmonTE/salmon/darwin/bin/salmon quant...' returned non-zero exit status 134.
```

**A.** You may have a problem to run `salmon` which is an essential tool for the pipeline. You may install [Threading Building Blocks library](https://www.threadingbuildingblocks.org/download) to solve the problem. If you are using [homebrew](https://brew.sh) then please use below command:

```
brew install tbb
```



## How to use it?

```
Usage:
    SalmonTE.py index [--ref_name=ref_name] (--input_fasta=fa_file) [--te_only]
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

## How to Cite?

```
@inbook{doi:10.1142/9789813235533_0016,
author = {Hyun-Hwan Jeong and Hari Krishna Yalamanchili and Caiwei Guo and Joshua M. Shulman and Zhandong Liu},
title = {An ultra-fast and scalable quantification pipeline for transposable elements from next generation sequencing data},
booktitle = {Biocomputing 2018},
chapter = {},
pages = {168-179},
doi = {10.1142/9789813235533_0016},
URL = {http://www.worldscientific.com/doi/abs/10.1142/9789813235533_0016},
eprint = {http://www.worldscientific.com/doi/pdf/10.1142/9789813235533_0016}
publisher = WORLD SCIENTIFIC
address = 
year = 2017
edition = 
}
```
