#!/usr/bin/env Rscript

suppressMessages(library(optparse))

option_list = list(
  make_option(c("-p", "--datapath"), action = "store", type = "character", default = NA,
              help = "Path to the directory containing necessary input data."),
  make_option(c("-o", "--resultpath"), action = "store", type = "character", default = NA,
              help = "Path to the directory for saving results.")
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
# dat_samples <- read.delim(file = file.path(path_data, "samples_QC.txt"), as.is = TRUE)
# samples <- sub("\\.txt$", "", dat_samples$File)
# n_samples <- nrow(dat_samples)

file_cnvr <- "cnvr_batch.txt"  ## with batch information
dt_cnvr_raw <- read.delim(file = file.path(path_data, file_cnvr), as.is = TRUE)

tbl_raw <- table(dt_cnvr_raw$chr, dt_cnvr_raw$batch)
dt_freq_raw <- as.data.frame(tbl_raw)
names(dt_freq_raw) <- c("chr", "batch", "Freq")
dt_freq_raw <- subset(dt_freq_raw, Freq != 0)

## initialize sample list using the information from the first CNVR
chr1   <- dt_freq_raw$chr[1]
batch1 <- dt_freq_raw$batch[1]
preds1 <- list.files(path = file.path(path_pred, paste0("chr_", chr1, "_batch_", batch1)),
                     pattern = ".rds")
dat1 <- readRDS( file = file.path(path_pred, paste0("chr_", chr1, "_batch_", batch1), preds1[1]) )

samples <- dat1$Sample_ID
n_samples <- length(samples)

cnvrs  <- c()

# row: CNVRs; column: samples
res_CN <- data.frame()
res_GQ <- data.frame()

for ( i in 1:nrow(dt_freq_raw) ) {
  
  chr1   <- dt_freq_raw$chr[i]
  batch1 <- dt_freq_raw$batch[i]

  preds1 <- list.files(path = file.path(path_pred, paste0("chr_", chr1, "_batch_", batch1)),
                       pattern = ".rds")
  
  cnvrs1 <- gsub("_pred.rds$", "", preds1, perl = T)
  cnvrs <- c(cnvrs, cnvrs1)
  
  res1_GQ <- matrix(nrow = length(cnvrs1), ncol = n_samples)
  rownames(res1_GQ) <- cnvrs1
  colnames(res1_GQ) <- samples
  res1_CN <- res1_GQ
  
  for (k in 1:length(preds1)) {
    pred1 <- preds1[k]
    cnvr1 <- cnvrs1[k]
    dat1 <- readRDS(file = file.path(path_pred, paste0("chr_", chr1, "_batch_", batch1), pred1))
    
    ## sort the results according to the order of samples  
    dat1 <- dat1[match(samples, dat1$Sample_ID), ]
    stopifnot( all(dat1$Sample_ID == samples) )
    
    res1_GQ[k, ] <- dat1$value_GQ
    res1_CN[k, ] <- dat1$CN_gatk_pred
  }
  
  res_GQ <- rbind(res_GQ, res1_GQ)
  res_CN <- rbind(res_CN, res1_CN)
}

stopifnot( all(rownames(res_GQ) == cnvrs) )
stopifnot( all(colnames(res_GQ) == samples) )
stopifnot( all(rownames(res_CN) == cnvrs) )
stopifnot( all(colnames(res_CN) == samples) )

mat_GQ <- as.matrix(res_GQ)
mat_CN <- as.matrix(res_CN)
rownames(mat_GQ) <- cnvrs
rownames(mat_CN) <- cnvrs
colnames(mat_GQ) <- samples
colnames(mat_CN) <- samples

## mark on successfully CNV-genotyped CNVRs
dt_cnvr_raw$genotype <- 0
dt_cnvr_raw$genotype[ dt_cnvr_raw$CNVR_ID %in% cnvrs ] <- 1

write.table(dt_cnvr_raw, 
            file = file.path(path_result, "cnvr_genotype.txt"),
            quote = F, row.names = F, sep = "\t")

write.table(data.frame(Sample_ID = samples), 
            file = file.path(path_result, "sample_genotype.txt"),
            quote = F, row.names = F, col.names = F, sep = "\t")

saveRDS(mat_GQ, file = file.path(path_result, "matrix_GQ.rds"))
saveRDS(mat_CN, file = file.path(path_result, "matrix_CN.rds"))

