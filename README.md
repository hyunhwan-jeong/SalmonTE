# SalmonTE

## Change Logs

* December 19, 2018: `SalmonTE` has been updated to `0.4` with some improvements!
  * Regarding [#14](https://github.com/hyunhwaj/SalmonTE/issues/14), now `SalmonTE` is coupled with the latest version of `snakemake`, so there is no more directory error. Furthermore, `SalmonTE` no longer runs with the older version of `snakemake`, please update the version of the `snakemake` package.
  * Regarding [#22](https://github.com/hyunhwaj/SalmonTE/issues/22), mappability for each file is now reported. (MAPPING_INFO.csv in the quantification output directory)
  * `test` function has been improved, and this will give better representations of the data.
  * Running of `SalmonTE` will not produce massive and enigmatic messages anymore.

* June 4, 2018: Now `SalmonTE` supports `fq` and `fq.gz` extensions.

* May 3, 2018: Update README.txt

* May 2, 2018: A bug fix regarding issue [#10](https://github.com/hyunhwaj/SalmonTE/issues/10)

* January 23, 2018: Added references of *Mus musculus*(mm) and *Danio rerio*(dr), added a function users to allow to build a customized index (`index` mode), and fixed a minor bug.

* December 25, 2017: Fixed a bug for single-end reads dataset.

* December 14, 2017: Fixed issue of the `macOS`. We have figured out there is a problem of old version of `snakemake`. If you already install the package, and the version is not above `4.0.0` (You can check it with `snakemake --version`) then please update version with below command:
`pip3 install snakemake --user --upgrade`

* November 27, 2017: source code of PSB manuscript is out now - our [manuscript](https://github.com/hyunhwaj/SalmonTE-manuscript/) and [response letter](https://github.com/hyunhwaj/SalmonTE-response).

* November 21, 2017: Version 0.2 is out, `SalmonTE` now supports two different type of expressions with `--exprtype` option. Use `--exprtype=TPM` if you want to use TPM values for the statistical analysis. If you want to run differential expression analysis, then I highly recommend to use `--exprtype=count`. [Here is the nice answer why](https://support.bioconductor.org/p/98820/).


## Notice

* June 5, 2018: Our recent paper about TE in Alzhimer's disease has been published in Cell Reports! Please read the paper to see how `SalmonTE` was succesfully applied! [(Cell Reports Link)](https://www.cell.com/action/showAbstract?pii=S2211-1247(18)30722-8)

## What is SalmonTE?
`SalmonTE` is an ultra-Fast and Scalable Quantification Pipeline of Transpose Element (TE) Abundances from Next Generation Sequencing Data. It comes with [Salmon](https://github.com/COMBINE-lab/salmon) which is a fast and accurate transcriptome quantification method. You can read the details of the pipeline and an example of real data study in [my recent published paper in PSB 2018](http://www.worldscientific.com/doi/10.1142/9789813235533_0016).

## What I need to run `SalmonTE`? Why I have to use it?
* You only need to have a set of FASTQ files and phenotype data. Furthermore, **`SalmonTE` automatically decided wether your dataset is paired-ends reads or not.** 
* conditions can be a numeric data or a categorical data. Based on the data type of the conditions of each sample, `SalmonTE` will run differential expression analysis or linear regression.
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
install.packages(c("tidyverse", "scales", "WriteXLS", "BiocManager"))
BiocManager::install("DESeq2", version = "3.8")
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
    SalmonTE.py test [--inpath=inpath] [--outpath=outpath] [--tabletype=tabletype] [--figtype=figtype] [--analysis_type=analysis_type] [--conditions=conditions]
    SalmonTE.py (-h | --help)
    SalmonTE.py --version

Options:
    -h --help     Show this screen.
    --version     Show version.
```

## An example of SalmonTE usage with command line 

### Running the `quant` mode to collect TE expressions

#### Parameters

* `--reference`: This will select a reference file, and should the species **identifier** of your data. We are currently supporting references of those species.
    * hs : *Homo Sapiens*
    * mm : *Mus musculus*
    * dm : *Drosophila melanogaster*
    * dr : *Danio rerio*
* `--outpath`: Path to generate quantification results. If you omit this, `SalmonTE_output` in the current path, and output will be stored in the path.
* `--exprtype`: The type of expression measure, and **TPM** or **count** are possible options. If you omit this, then "TPM" is the input of the parameter.
* `--num_threads`: This has to be an integer, and this parameter will configure how many threads will use for the parallization.

After you put your parameters, you can put the directory which includes a list of **FASTQ** files, 

```
SalmonTE.py quant --reference=hs example
```

Or, you can put the list of files like below.

```
SalmonTE.py quant --reference=hs example/CTRL_1_R1.fastq.gz example/CTRL_2_R1.fastq.gz          
```

### Running `test` mode to perform statistical test

Before you run test mode, you should modify `control.csv` file which is stored in the `outpath`. Here are examples of the proper modifications:

For the differential expression analysis, change the file as below. 
**Important**: The control samples has to be labeled as `control`. Other labels will cause errors.

```
SampleID,control
FASTQ1,control
FASTQ2,control
FASTQ3,treatment
FASTQ4,treatment
```

For the regression analysis, 

```
SampleID,control
FASTQ1,1.5
FASTQ2,2.1
FASTQ3,3.8
FASTQ4,9.5
```

Once the conditions of every sample has been filled, we can run the test mode like the example commnad-line below:

* `--inpath`: This should be the path which contains output of `quant` mode. 
* `--outpath`: This will be the path to store all outputs for the mode.
* `--tabletype`: The file format of the tables, `csv`, `tsv`, and `xls` are supported. If you omit this, then `xls` formatted file will be generated.
* `--tabletype`: The file format of the figures, `png`, `tiff`, `jpg`, and `pdf` are supported. If you omit this, then `pdf` formated files will be generated.
* `--analysis_type`: The type of the analysis, and **DE** (for a differential analysis) or **LM** (for a linear regression analysis) are possible options. If you omit this, then "DE" is the input of the parameter.
* `--conditions`: The list of conditions will be considered when **DE** has been selected for `--analysis_type.` The input needs to contain two different conditions (written in non-white space characters) and each condition are separated by `,`, and **no white-space characters** are not allowed for the input. i.e. `SalmonTE` does not care input such as `control , treatment`. The first condition of the input will be considered as a normal condition (e.g. healthy condition, wild-type mice) in the study, and the later will be considered as another condition which you are interested (e.g. knock-out mice, treatment).

```
SalmonTE.py test --inpath=SalmonTE_output --outpath=SalmonTE_statistical_test --tabletype=csv --figtype=png --analysis_type=DE --conditions=control,treatment
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
