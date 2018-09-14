#!/usr/bin/env Rscript

suppressMessages( library(optparse) )
# before running this script 
# you need to regenotype CNVRs in dat_cnvrs_regenotype.rds
# which are generated from step4.clean.results

option_list <- list(
  make_option(c("-b", "--matrixbeforerefine"), action = "store", type = "character", default = NA,
              help = "matrix of copy number before refinement."),
  make_option(c("-f", "--matrixrefine"), action = "store", type = "character", default = NA, 
              help = "matrix from refinement CNVR."),
  make_option(c("-p", "--refinepath"), action = "store", type = "character", default = NA,
              help = "path saving cnvrs_keep_after_refine.rds."),
  make_option(c("-o", "--output"), action = "store", type = "character", default = NA,
              help = "path to saving final matrix")
)

opt <- parse_args(OptionParser(option_list = option_list))
pars <- c(opt$matrixbeforerefine, opt$matrixrefine, 
          opt$refinepath, opt$output)
if ( any(is.na(pars)) ) {
  stop("All parameters must be supplied. (--help for detail)")
}

path_matrix_before_refine <- opt$matrixbeforerefine
path_matrix_refine <- opt$matrixrefine
path_refine <- opt$refinepath
path_output <- opt$output


mat_CN_before_refine <- readRDS( file = file.path(path_matrix_before_refine, "mat_CN.rds"))
mat_GQ_before_refine <- readRDS( file = file.path(path_matrix_before_refine, "mat_GQ.rds"))

cnvrs <- rownames( mat_CN_before_refine )
samples <- colnames( mat_CN_before_refine )

# keep cnvrs after refinement
cnvrs_keep_final <- readRDS( file = file.path(path_refine, "cnvrs_keep_after_refine.rds"))

mat_CN_keep <- mat_CN_before_refine[cnvrs_keep_final, ]
mat_GQ_keep <- mat_GQ_before_refine[cnvrs_keep_final, ]

# regenotype refinement CNVR ----------------------------------------------

mat_CN_refine <- readRDS( file = file.path(path_matrix_refine, "mat_CN.rds"))
mat_GQ_refine <- readRDS( file = file.path(path_matrix_refine, "mat_GQ.rds"))

samples_refine <- colnames( mat_CN_refine )
stopifnot( sum(samples_refine %in% samples) == length(samples))

mat_CN_refine <- mat_CN_refine[, samples]
mat_GQ_refine <- mat_GQ_refine[, samples]


mat_CN_final <- rbind( mat_CN_keep, mat_CN_refine )
mat_GQ_final <- rbind( mat_GQ_keep, mat_GQ_refine )

saveRDS( mat_CN_final, file = file.path(path_output, "mat_CN_final.rds"))
saveRDS( mat_GQ_final, file = file.path(path_output, "mat_GQ_final.rds"))



