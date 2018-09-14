#!/usr/bin/env Rscript

suppressMessages(library(optparse))

option_list <- list(
  make_option(c("-i", "--datapath"), action = "store", type = "character", default = NA,
              help = "path contain all running needed data"),
  make_option(c("-o", "--resultpath"), action = "store", type = "character", default = NA,
              help = "path save all results."),
  make_option(c("-c", "--cutofffreq"), action = "store", type = "double", default = NA,
              help = "cutoff of CNVR frequence to select refinement CNVRs.")
)


opt = parse_args(OptionParser(option_list = option_list))
pars = c(opt$datapath, opt$resultpath, opt$cutofffreq)

if ( any(is.na(pars)) ) {
  stop("All three parameters must be supplied. (--help for detail)")
}

cutoff_freq <- as.numeric( opt$cutofffreq)
path_data   <- opt$datapath
path_result <- opt$resultpath

path_output <- file.path(path_result, "cnvr_refinement")
if (!dir.exists(paths = path_output) ) dir.create(path = path_output, showWarnings = F, recursive = T)

# here using the matrix of copy number from regenotype step
mat_regenotype <- readRDS( file = file.path(path_result, "mat_CN.rds"))
n.sample <- ncol( mat_regenotype )
n.CNVR   <- nrow( mat_regenotype )

cnvrs <- rownames( mat_regenotype )

freqs_CNVR <- unlist( lapply(1:n.CNVR, FUN = function(i) {
  v1 <- as.integer( mat_regenotype[i, ])
  n1 <- sum( v1 %in% c(0, 1, 3))
  n1
}))

idxs.refine <- which( freqs_CNVR >= n.sample*cutoff_freq)

dat_freq <- data.frame(CNVR_ID = cnvrs, Freq = freqs_CNVR,
                       stringsAsFactors = F)

cnvrs_refine <- cnvrs[ idxs.refine ]
cnvrs_keep   <- cnvrs[ -idxs.refine ]

saveRDS( cnvrs_refine, file = file.path(path_output, "cnvrs_refine.rds"))
saveRDS( cnvrs_keep, file = file.path(path_output, "cnvrs_keep.rds"))

file_cnvr  <- "cnvr_batch.txt"  ## with batch information
dat_cnvrs  <- read.delim(file = file.path(path_data, file_cnvr), as.is = TRUE)
dat_cnvrs  <- merge( dat_cnvrs, dat_freq)
stopifnot( nrow(dat_cnvrs) == nrow(dat_freq) )

dat_cnvrs_refine <- subset( dat_cnvrs, CNVR_ID %in% cnvrs_refine )
dat_cnvrs_keep   <- subset( dat_cnvrs, CNVR_ID %in% cnvrs_keep )

saveRDS( dat_cnvrs_keep, file = file.path(path_output, "dat_cnvrs_keep.rds"))
saveRDS( dat_cnvrs_refine, file = file.path(path_output, "dat_cnvrs_refine.rds"))
