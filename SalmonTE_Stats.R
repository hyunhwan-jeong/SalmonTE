setwd("~/Sandbox/SalmonTE/")
library(getopt)
library(tidyverse)
library(DESeq2)
library(tximport)
library(cowplot)
library(scales)
library(WriteXLS)

write.result <- function(dat) {
  res <- dat$res
  sheet.fmt <- dat$sheet.fmt
  path <- dat$path
  if(sheet.fmt == "tsv") {
    write.table(df.res, file=paste0(path,"/results.tsv"), sep="\t")
  } else if (sheet.fmt == "csv") {
    write.csv(df.res, file=paste0(path,"/results.csv"))
  } else {
    WriteXLS(df.res, ExcelFileName=paste0(path,"/results.xls"), row.names=T)
  }
}


do.deseq2 <- function(dat) {
  tpm <- dat$tpm
  phenoData <- dat$phenoData
  dds <- DESeqDataSetFromMatrix(countData = ceiling(tpm),
                                colData = phenoData,
                                design = ~ phenotype)
  
  dds <- dds[ rowSums(counts(dds)) > 1, ]
  dds <- DESeq(dds)
  res <- results(dds)
  data.frame(res)
}

do.lm <- function(dat) {
  tpm <- dat$tpm
  phenoData <- dat$phenoData
  tidy.tpm <- tpm %>% 
    cbind(name = rownames(tpm)) %>%
    gather(sample, tpm, -name) %>% left_join(phenoData %>% mutate(name = rownames(phenoData) ), by = c("sample"="name"))
  
  # TIP: values in a column must be atomic, can't have a vector
  res <- tidy.tpm %>% group_by(name) %>% 
    summarise(tpm = list(tpm), pheno = list(phenotype)) %>%
    group_by(name) %>%
    mutate( lm = list(summary(lm(unlist(pheno)~unlist(tpm)))),
            baseMean = mean(unlist(tpm))) %>%
    mutate( p.value = tryCatch({lm[[1]]$coefficients[2,4]}, error = function(e) NA),
            b.value = tryCatch({lm[[1]]$coefficients[2,3]}, error = function(e) NA)) %>%
    select(name, baseMean, b.value, p.value )
  res$padj <- p.adjust(res$p.value, method="fdr")
  res
}

SalmonTE <- function(tpm, phenoData, sheet.fmt = "csv", path = ".") {
  dat <- list(
    tpm = tpm, phenoData = phenoData
  )
  if( "control" %in% phenoData$phenotype || "case" %in% phenoData$phenotype ) {
    phenoData$phenotype <- factor(phenoData$phenotype)
    dat$res <- do.deseq2(dat)
    dat$y.name <- "log2FoldChange"
  }
  else {
    dat$res <- do.lm(dat)
    dat$y.name <- "b.value"
  }
  dat$sheet.fmt <- sheet.fmt
  dat$path <- path
  dat
}

# a code to draw MA-plot
draw.MAplot <- function(dat) {
  dat$res %>%
    filter( !is.na(padj) ) %>%
    mutate( sig = ifelse(padj < 0.05, "DE", "NC" )) %>%
    ggplot( aes_string(x="baseMean", y=dat$y.name))+
    geom_hline(yintercept = 0, col = "red", alpha=0.5) +
    geom_point( aes(colour=sig) ) + 
    scale_colour_manual(values = c("red", "black"), limits = c("DE", "NC")) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)) ) +
    theme(legend.position = "none")
}

tpm <- read.csv("SalmonTE_output/TPM.csv", row.names="TE")
phenoData <- read.csv("SalmonTE_output/phenotype.csv", 
                      row.names = "SampleID")

dat <- SalmonTE(tpm, phenoData)
draw.MAplot(dat)
write.result(dat)
