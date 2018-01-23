library(tidyverse)
library(DESeq2)
library(tximport)
library(cowplot)
library(scales)
library(WriteXLS)

write.results <- function(dat) {
  res <- dat$res
  summary <- dat$summary
  sheet.fmt <- dat$sheet.fmt
  path <- dat$path
  if(sheet.fmt == "tsv") {
    write.table(res, file=file.path(path,"results.tsv"), sep="\t",  row.names = F)
    write.table(summary$clade$table, file=file.path(path,"results.clade.tsv"), sep="\t", row.names = F)
    write.table(summary$class$table, file=file.path(path,"results.class.tsv"), sep="\t", row.names = F)
  } else if (sheet.fmt == "csv") {
    write.csv(res, file=file.path(path,"results.csv"), row.names = F)
    write.csv(summary$clade$table, file=file.path(path,"results.clade.csv"), row.names = F)
    write.csv(summary$class$table, file=file.path(path,"results.class.csv"), row.names = F)
  } else {
    WriteXLS(res, ExcelFileName=file.path(path,"results.xls"), row.names=F)
    WriteXLS(summary$clade$table, ExcelFileName=file.path(path,"results.clade.xls"), row.names=F)
    WriteXLS(summary$class$table, ExcelFileName=file.path(path,"results.class.xls"), row.names=F)
  }
}

write.figures <- function(dat) {
  path <- dat$path
  ext <- dat$fig.fmt
  save_plot(file.path(path, paste0("ma.plot.", ext)), dat$summary$ma.plot)
  save_plot(file.path(path, paste0("clade.", ext)), dat$summary$clade$figure)
  save_plot(file.path(path, paste0("class.", ext)), dat$summary$class$figure)
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

do.summary <- function(dat) {
  calc.stat <- function(group) {
    tmp <- data.frame(group=dat$res[,group], 
                      value=dat$res[,dat$y.name]) 
    colnames(tmp) <- c("group", "value")
    tmp <- tmp %>%
      filter(!is.na(value))
    tb <- tmp %>%
      group_by(group) %>% 
      summarise( n = n(), 
                 mean = mean(value, na.rm=T), 
                 sd = sd(value, na.rm=T),
                 p.value = tryCatch({t.test(value)$p.value}, 
                                    error = function(e) NA ) )
    fig <- tmp %>% 
      ggplot(aes(x=group, y=value)) +
      geom_boxplot(aes(fill=group)) +
      geom_jitter( width = 0.5) + theme(legend.position = "none") +
      xlab(group) + ylab(dat$y.name)
    return(list(table = tb, figure = fig))
  }
  dat$summary <- list(class = calc.stat("class"),
                      clade = calc.stat("clade"))
  return(dat)
}

# a code to draw MA-plot
draw.MAplot <- function(dat) {
  return(dat$res %>%
    filter( !is.na(padj) ) %>%
    mutate( sig = ifelse(padj < 0.05, "DE", "NC" )) %>%
    ggplot( aes_string(x="baseMean", y=dat$y.name))+
    geom_hline(yintercept = 0, col = "red", alpha=0.5) +
    geom_point( aes(colour=sig) ) + 
    scale_colour_manual(values = c("red", "black"), limits = c("DE", "NC")) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)) ) +
    theme(legend.position = "none"))
}

SalmonTE <- function(tpm, phenoData, annotation, sheet.fmt = "csv", fig.fmt = "pdf", path = ".") {
  dat <- list(
    tpm = tpm, phenoData = phenoData
  )
  if( "control" %in% phenoData$phenotype || "case" %in% phenoData$phenotype ) {
    phenoData$phenotype <- factor(phenoData$phenotype)
    dat$res <- do.deseq2(dat)
    dat$res <- cbind(data.frame(name=rownames(dat$res)), dat$res)
    dat$y.name <- "log2FoldChange"
  } else {
    dat$res <- do.lm(dat)
    dat$y.name <- "b.value"
  }
  dat$annotation <- annotation
  dat$sheet.fmt <- sheet.fmt
  dat$path <- path
  dat$fig.fmt <- fig.fmt
  dat$res <- left_join(dat$res, dat$annotation, by="name")
  nc <- ncol(dat$res)
  dat$res <- dat$res[, c(1, nc-1, nc, which(!(1:nc %in% c(1,nc-1,nc))))]

  dat <- do.summary(dat)
  dat$summary$ma.plot <- draw.MAplot(dat)
  dat
}

GenerateOutput <- function(dat) {
  write.results(dat)
  write.figures(dat)
  save(dat, file = file.path(dat$path, "data.Rdata"))
}

args <- commandArgs(T)

tpm <- read.csv(file.path(args[1], "EXPR.csv"), row.names="TE")
phenoData <- read.csv(file.path(args[1], "phenotype.csv"), 
                      row.names = "SampleID")
annotation <- read.csv(file.path(args[1], "clades.csv"))
dat <- SalmonTE(tpm, phenoData, annotation, args[2], args[3], args[4])
GenerateOutput(dat)
