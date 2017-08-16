setwd("~/Sandbox/SalmonTE/")
library(getopt)
library(tidyverse)
library(DESeq2)
library(tximport)
library(cowplot)
library(scales)
library(WriteXLS)
tpm <- read.csv("SalmonTE_output/TPM.csv", row.names="TE")
phenoData <- read.csv("SalmonTE_output/phenotype.csv", 
                      row.names = "SampleID")

dds <- DESeqDataSetFromMatrix(countData = ceiling(tpm),
                              colData = phenoData,
                              design = ~ phenotype)

dds <- dds[ rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)
res <- results(dds)
df.res <- data.frame(res)
summary(res)

# declations to write results
write.table(res, file="SalmonTE_output/results.tsv", sep="\t")
write.csv(res, file="SalmonTE_output/results.csv")
WriteXLS(df.res, ExcelFileName="SalmonTE_output/results.xls", row.names=T)

# a code to draw MA-plot
df.res %>% 
  mutate( sig = ifelse(padj < 0.05, "DE", "NC" )) %>%
  ggplot( aes(x=baseMean, y=log2FoldChange))+
  geom_hline(yintercept = 0, col = "red", alpha=0.5) +
  geom_point( aes(colour=sig) ) + 
  scale_colour_manual(values = c("red", "black"), limits = c("DE", "NC")) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)) ) +
  theme(legend.position = "none")
