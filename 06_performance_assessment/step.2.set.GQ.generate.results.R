#!/usr/bin/env Rscript

suppressMessages(library(optparse))
suppressMessages(library(plyr))

option_list <- list(
  make_option(c("-n", "--matrixCN"), action = "store", default = NA,type = "character",
              help = "Path to matrix of copy number (CN)"),
  make_option(c("-g", "--matrixGQ"), action = "store", default = NA,type = "character",
              help = "Path to matrix of genotyping quality (GQ) score."),
  make_option(c("-c", "--cnvrfile"), action = "store", default = NA, type = "character",
              help = "Path to CNVR information after boundary refinement."),
  make_option(c("-o", "--resultpath"), action = "store", default = NA,type = "character",
              help = "Path to directory for saving assessment results."),
  make_option(c("-s", "--gqscore"), action = "store", default = NA, type = "integer",
              help = "Set GQ score threshold.")
)

opt <- parse_args(OptionParser(option_list = option_list))
pars <- c(opt$matrixCN, opt$matrixGQ, opt$cnvrfile,
          opt$resultpath, opt$gqscore)

if (any(is.na(pars))) {
  stop("All required parameters must be supplied. (--help for detail)")
}

file_matrixcn <- opt$matrixCN
file_matrixgq <- opt$matrixGQ
file_cnvr     <- opt$cnvrfile
path_result   <- opt$resultpath
gqscore       <- as.numeric(opt$gqscore)

matrix_CN <- readRDS(file = file_matrixcn)
matrix_gq <- readRDS(file = file_matrixgq)
dat_cnvr  <- read.delim(file = file_cnvr, check.names = FALSE, as.is = TRUE)

# main --------------------------------------------------------------------

idxs.nocall = which(matrix_gq < gqscore)

if (length(idxs.nocall) >= 1) matrix_CN[idxs.nocall] = -9

cnvrs   <- rownames( matrix_CN )
samples <- colnames( matrix_CN )

n_cnvr   <- nrow(matrix_CN)
n_sample <- ncol(matrix_CN)

## cnvr freq
list_freqs_cnvr <- lapply(1:n_cnvr, FUN = function(k) {
  v1 <- as.vector(matrix_CN[k, ])
  data.frame(n = length(v1), 
             n0 = sum(v1 == 0),
             n1 = sum(v1 == 1),
             n2 = sum(v1 == 2),
             n3 = sum(v1 == 3),
             n_nocall = sum(v1 == -9))
})

freqs_cnvr <- do.call(rbind, list_freqs_cnvr)

dat_freqs_cnvr <- data.frame(freqs_cnvr, stringsAsFactors = F, check.names = F)
dat_freqs_cnvr$CNVR_ID <- cnvrs

dat_freqs_cnvr$callRate <- (dat_freqs_cnvr$n0 + dat_freqs_cnvr$n1 + dat_freqs_cnvr$n2 + dat_freqs_cnvr$n3)/dat_freqs_cnvr$n
dat_freqs_cnvr$freq     <- (dat_freqs_cnvr$n0 + dat_freqs_cnvr$n1 + dat_freqs_cnvr$n3)/dat_freqs_cnvr$n

idxs_cnvr_filter <- which(dat_freqs_cnvr$freq == 0)

if (length(idxs_cnvr_filter) >= 1) {
  dat_freqs_cnvr <- dat_freqs_cnvr[-idxs_cnvr_filter, ]
}
dat_cnvr_final <- merge(dat_cnvr, dat_freqs_cnvr)
dat_cnvr_final <- dat_cnvr_final[order(dat_cnvr_final$chr, 
                                       dat_cnvr_final$arm, 
                                       dat_cnvr_final$posStart,
                                       dat_cnvr_final$posEnd), ]

dat_cnvr_final <- dat_cnvr_final[, c("CNVR_ID", "chr", "arm", "posStart", "posEnd", "start_snp", "end_snp", 
                                     "n", "n0", "n1", "n2", "n3", "n_nocall", "callRate", "freq")]
dat_cnvr_final <- rename(dat_cnvr_final, c("start_snp"="snpStart", "end_snp"="snpEnd"))

cat(nrow(dat_cnvr_final), "CNVRs remains from", nrow(n_cnvr), "CNVRs after GQ cut-off.\n")

write.table(dat_cnvr_final, file = file.path(path_result, "cnvr_after_GQ.txt"),
            sep = "\t", row.names = F, col.names = T, quote = F)

matrix_CN_final <- matrix_CN[dat_cnvr_final$CNVR_ID, ]
saveRDS(matrix_CN_final, file = file.path(path_result, "matrix_CN_after_GQ.rds"))

# sample information ------------------------------------------------------

list_samples_info <- lapply(1:ncol(matrix_CN_final), FUN = function(k) {
  v1 <- as.vector(matrix_CN_final[, k])
  data.frame(n = length(v1),
             n0 = sum(v1 == 0),
             n1 = sum(v1 == 1),
             n2 = sum(v1 == 2),
             n3 = sum(v1 == 3),
             n_nocall = sum(v1 == -9))
})

samples_info <- do.call(rbind, list_samples_info)
samples_info <- data.frame(samples_info, stringsAsFactors = F, check.names = F)
samples_info$Sample_ID <- samples

samples_info$callRate <- (samples_info$n0 + samples_info$n1 + samples_info$n2 + samples_info$n3)/samples_info$n
samples_info$freq     <- (samples_info$n0 + samples_info$n1 + samples_info$n3)/samples_info$n
samples_info <- samples_info[, c("Sample_ID", "callRate", "freq", "n", "n0", "n1", "n2", "n3", "n_nocall")]

write.table(samples_info, file = file.path(path_result, "sample_after_GQ.txt"),
            sep = "\t", row.names = F, col.names = T, quote = F)

