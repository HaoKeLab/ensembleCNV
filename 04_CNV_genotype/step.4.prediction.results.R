#!/usr/bin/env Rscript

suppressMessages(library(optparse))

option_list = list(
  make_option(c("-p", "--datapath"), action = "store", type = "character", default = NA,
              help = "path contain all running needed data"),
  make_option(c("-o", "--resultpath"), action = "store", type = "character", default = NA,
              help = "path save all results")
)

opt = parse_args(OptionParser(option_list = option_list))
pars = c(opt$datapath, opt$resultpath)

if ( any(is.na(pars)) ) {
  stop("All parameters must be supplied. (--help for detail)")
}

path_data   <- opt$datapath
path_result <- opt$resultpath

path_pred <- file.path(path_result, "pred")

# number of samples
dat_samples <- read.delim(file = file.path(path_data, "samples_QC.txt"), as.is = TRUE)
n_samples <- nrow(dat_samples)

file_cnvr <- "cnvr_batch.txt"  ## with batch information
dt_cnvr_raw <- read.delim(file = file.path(path_data, file_cnvr), as.is = TRUE)
tbl_raw <- table(dt_cnvr_raw$chr, dt_cnvr_raw$batch)
dt_freq_raw <- as.data.frame(tbl_raw)
names(dt_freq_raw) <- c("chr", "batch", "Freq")

dt_freq_raw <- subset(dt_freq_raw, Freq != 0)

res_CN  <- data.frame()
res_GQ  <- data.frame()
cnvrs   <- c()
samples <- c()
flag_samples <- TRUE

# rowname = "CNVR_ID"; colname = "sample_ID"
for ( i in 1:nrow(dt_cnvr_raw) ) {
  
  chr1   <- dt_cnvr_raw$chr[i]
  batch1 <- dt_cnvr_raw$batch[i]

  preds1 <- list.files(path = file.path(path_pred, paste0("chr_", chr1, "_batch_", batch1)),
                       pattern = ".rds")
  cnvrs1 <- gsub("_pred.rds$", "", preds1, perl = T)
  
  cnvrs <- c(cnvrs, cnvrs1)
  res1_GQ <- matrix(nrow = length(cnvrs), ncol = n_samples)
  res1_CN <- matrix(nrow = length(cnvrs), ncol = n_samples)
  for (k in 1:length(preds1)) {
    pred1 <- preds1[k]
    cnvr1 <- cnvrs1[k]
    dat1 <- readRDS(file = file.path(path_pred, paste0("chr_", chr1, "_batch_", batch1), pred1))
    
    if ( flag_samples) samples <- dat1$Sample_ID
    
    dat1 <- dat1[match(samples, dat1$Sample_ID), ]
    stopifnot( all(dat1$Sample_ID == samples) )
    
    res1_GQ[k, ] <- dat1$value_GQ
    res1_CN[k, ] <- dat1$CN_gatk_pred
  }
  
  res_GQ <- rbind(res_GQ, res1_GQ)
  res_CN <- rbind(res_CN, res1_CN)
}

mat_GQ <- as.matrix(res_GQ)
mat_CN <- as.matrix(res_CN)

rownames(mat_GQ) <- cnvrs
rownames(mat_CN) <- cnvrs

colnames(mat_GQ) <- samples
colnames(mat_CN) <- cnvrs

saveRDS(cnvrs, file = file.path(path_result, "cnvrs_all_predict.rds"))
saveRDS(samples, file = file.path(path_result, "samples_all_predict.rds"))
saveRDS(mat_GQ, file = file.path(path_result, "matrix_GQ.rds"))
saveRDS(mat_CN, file = file.path(path_result, "matrix_CN.rds"))

