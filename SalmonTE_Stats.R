message("Step 1: Loading required libraries...")
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(WriteXLS))
suppressPackageStartupMessages(library(ggforce))
suppressPackageStartupMessages(library(ggsci))

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
  ggsave(file.path(path, paste0("ma.plot.", ext)), dat$summary$ma.plot)
  ggsave(file.path(path, paste0("clade.", ext)), dat$summary$clade$figure)
  ggsave(file.path(path, paste0("class.", ext)), dat$summary$class$figure)
}

do.deseq2 <- function(dat) {
  count <- dat$count
  col_data <- dat$col_data
  dds <- DESeqDataSetFromMatrix(countData = round(count),
                                colData = col_data,
                                design = ~condition)
  dds <- dds[ rowSums(counts(dds)) > 1, ]
  dds <- DESeq(dds)
  res <- results(dds)
  df_res <- data.frame(res)
  dat$res <- df_res %>% rownames_to_column("name")
  dat$dds <- dds
  
  res$padj < 0.01 & abs(res$log2FoldChange) > 0.5
  dat
}


do.lm <- function(dat) {
  count <- dat$count
  col_data <- dat$col_data
  tidy.count <- count %>% 
    cbind(name = rownames(count)) %>%
    gather(sample, count, -name) %>% 
      left_join(col_data %>% rownames_to_column("sample"), by = "sample")    
  # TIP: values in a column must be atomic, can't have a vector
  res <- tidy.count %>% group_by(name) %>% 
    summarise(count = list(count), 
              condition = list(condition)) %>%
    group_by(name) %>%
    mutate( lm = list(summary(lm(unlist(condition)~unlist(count)))),
            baseMean = mean(unlist(count))) %>%
    mutate( p.value = tryCatch({lm[[1]]$coefficients[2,4]}, error = function(e) NA),
            b.value = tryCatch({lm[[1]]$coefficients[2,3]}, error = function(e) NA)) %>%
    select(name, baseMean, b.value, p.value)
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
      geom_boxplot(width=0.5) + 
      geom_jitter(width=0.1, alpha=0.5) + theme(legend.position = "none") +
      geom_hline(yintercept = 0, color='red') +
      xlab(group) + ylab(dat$y.name) + scale_fill_npg() + theme_minimal() + 
      theme(text = element_text(size = 18)) +
      theme(legend.position = "none")
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
    theme(legend.position = "none", text = element_text(size = 18))  + theme_minimal())
}

 SalmonTE <- function(count, col_data, annotation,
                      analysis,
                      condition_level,
                      sheet.fmt = "csv", fig.fmt = "pdf", path = ".") {
  dat <- list(
    count = count, 
    col_data = col_data
  )

  if( analysis == "DE" ) {
    if(!is.null(condition_level)) {
      dat$col_data$condition <- factor(dat$col_data$condition, level = condition_level)
    }
    dat <- do.deseq2(dat)
    dat$y.name <- "log2FoldChange"
  } else {
    dat$res <- do.lm(dat)
    dat$y.name <- "b.value"
  }
  dat$annotation <- annotation
  dat$sheet.fmt <- sheet.fmt
  dat$path <- path
  dat$fig.fmt <- fig.fmt
  dat$res <- suppressWarnings(left_join(dat$res, dat$annotation, by="name"))
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

message("Step 2: Loading input data...")
args <- commandArgs(T)

analysis <- args[5]
condition_level <- NULL
if(analysis == "DE" && !is.na(args[6])) {
  condition_level <- str_split(args[6], ",", simplify = T)
}
count <- read.csv(file.path(args[1], "EXPR.csv"), row.names="TE", check.names = FALSE)
col_data <- read.csv(file.path(args[1], "condition.csv"), row.names = "SampleID")
annotation <- read.csv(file.path(args[1], "clades.csv"))
message(sprintf("Step 3: Running the %s analysis...", analysis))
dat <- SalmonTE(count, col_data, annotation, analysis, condition_level, args[2], args[3], args[4])

message(sprintf("Step 4: Generating output...", analysis))
suppressMessages(GenerateOutput(dat))
message(sprintf("Step 5: The statistical analysis has been completed. Please check '%s' directory to see the analysis result!", args[4]))
