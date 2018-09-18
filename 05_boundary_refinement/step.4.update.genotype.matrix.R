#!/usr/bin/env Rscript

suppressMessages( library(optparse) )
# before running this script 
# you need to regenotype CNVRs in cnvr_regenotype_after_refine.txt
# which are generated from step3.clean.results.R

option_list <- list(
  make_option(c("-b", "--matrixbeforerefine"), action = "store", type = "character", default = NA,
              help = "Path to CN and GQ matrices generated in first round of 
                      CNV genotyping step before boundary refinement."),
  make_option(c("-f", "--matrixrefine"), action = "store", type = "character", default = NA, 
              help = "Path to CN and GQ matrices generated in CNV regenotyping for 
                      CNVRs with updated boundaries after refinement as well as CNVR information."),
  make_option(c("-p", "--refinepath"), action = "store", type = "character", default = NA,
              help = "Path to cnvr_kept_after_refine.txt."),
  make_option(c("-o", "--output"), action = "store", type = "character", default = NA,
              help = "Path to the directory for saving final CN and GQ matrices")
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

mat_CN_before_refine <- readRDS( file = file.path(path_matrix_before_refine, "matrix_CN.rds"))
mat_GQ_before_refine <- readRDS( file = file.path(path_matrix_before_refine, "matrix_GQ.rds"))

cnvrs   <- rownames( mat_CN_before_refine )
samples <- colnames( mat_CN_before_refine )

# keep cnvrs after refinement
dat_cnvr_keep <- read.delim( file = file.path(path_refine, "cnvr_kept_after_refine.txt"), as.is = TRUE)

mat_CN_keep <- mat_CN_before_refine[dat_cnvr_keep$CNVR_ID, ]
mat_GQ_keep <- mat_GQ_before_refine[dat_cnvr_keep$CNVR_ID, ]

# regenotyped CNVRs with updated boundaries ----------------------------------------------
dat_cnvr_refine <- read.delim( file = file.path(path_matrix_refine, "cnvr_genotype.txt"), as.is = TRUE)

mat_CN_refine <- readRDS( file = file.path(path_matrix_refine, "matrix_CN.rds"))
mat_GQ_refine <- readRDS( file = file.path(path_matrix_refine, "matrix_GQ.rds"))

samples_refine <- colnames( mat_CN_refine )
stopifnot( sum(samples_refine %in% samples) == length(samples))

mat_CN_refine <- mat_CN_refine[, samples]
mat_GQ_refine <- mat_GQ_refine[, samples]

## final results
mat_CN_final <- rbind( mat_CN_keep, mat_CN_refine )
mat_GQ_final <- rbind( mat_GQ_keep, mat_GQ_refine )

common.cols <- intersect(names(dat_cnvr_keep), names(dat_cnvr_refine))
common.cols <- setdiff(common.cols, c("batch", "genotype", "identicalID"))
dat_cnvr <- rbind(dat_cnvr_keep[, common.cols], dat_cnvr_refine[, common.cols])

saveRDS( mat_CN_final, file = file.path(path_output, "matrix_CN_final.rds"))
saveRDS( mat_GQ_final, file = file.path(path_output, "matrix_GQ_final.rds"))

write.table(dat_cnvr,
            file = file.path(path_output, "cnvr_final.txt"),
            quote = F, row.names = F, sep = "\t")

