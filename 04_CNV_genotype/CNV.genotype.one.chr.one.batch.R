#!/usr/bin/env Rscript --vanilla

# load packages
suppressMessages({
  require(optparse)
  require(dplyr)
  require(mixtools)
  require(ggplot2)
  require(cowplot)
  require(plyr)
  require(modeest)
  require(mclust)
  require(gridExtra)
  require(pheatmap)
  require(RColorBrewer)
})

option_list = list(
  make_option(c("-c", "--chr"), action = "store", type = "character", default = NA,
              help = "which chromosome."),
  make_option(c("-b", "--batch"), action = "store", type = "character", default = NA,
              help = "which batch."),
  make_option(c("-t", "--type"), action = "store", type = "character", default = NA,
              help = "Job submission type (0 - initial submission, 1 - resubmission of failed jobs)"),
  make_option(c("-p", "--datapath"), action = "store", type = "character", default = NA,
              help = "path contain all running needed data"),
  make_option(c("-o", "--resultpath"), action = "store", type = "character", default = NA,
              help = "path save all results"),
  make_option(c("-m", "--matrixpath"), action = "store", type = "character", default = NA,
              help = "path save all LRR and BAF chromosome based data"),
  make_option(c("-s", "--sourcefile"), action = "store", type = "character", default = NA,
              help = "path contain all scripts need to be soucred")
)

opt = parse_args(OptionParser(option_list = option_list))
pars = c(opt$chr, opt$batch, opt$type, opt$datapath, opt$resultpath, opt$matrixpath, opt$sourcefile)

if ( any(is.na(pars)) ) {
  stop("All three parameters must be supplied. (--help for detail)")
}

# opt <- list()
# opt$chr <- "1"
# opt$batch <- "1"
# opt$type <- "0"
# opt$datapath <- "/sc/orga/projects/haok01a/chengh04/paper/ensembleCNV_code_test/test/04_CNV_genotype/data/"
# opt$resultpath <- "/sc/orga/projects/haok01a/chengh04/paper/ensembleCNV_code_test/test/04_CNV_genotype/mid_res/"
# opt$matrixpath <- "/sc/orga/projects/haok01a/chengh04/paper/ensembleCNV_code_test/test/01_initial_call/finalreport_to_matrix_LRR_and_BAF/RDS/"
# opt$sourcefile <- "/sc/orga/projects/haok01a/chengh04/paper/ensembleCNV_code_test/code/04_CNV_genotype/scripts/"

chr1   <- as.integer( opt$chr )
batch1 <- as.integer( opt$batch )
type1  <- as.integer( opt$type )

path_data   <- opt$datapath
path_result <- opt$resultpath
path_matrix <- opt$matrixpath
path_sourcefile <- opt$sourcefile

if ( type1 != 1 & type1 != 0) {
  stop("Job submission type must be 0 or 1. (--help for detail)")
}
## print out parameters
cat("Processing chr:", chr1, "batch:", batch1, "type:", type1, "\n")

# source all the functions used in the pipeline
source(file = file.path(path_sourcefile, 'fun_BAF.R'))
source(file = file.path(path_sourcefile, 'fun_gatk.R'))
source(file = file.path(path_sourcefile, 'fun_LRR.R'))
source(file = file.path(path_sourcefile, 'fun_models.R')) 
source(file = file.path(path_sourcefile, 'fun_plot_steps.R'))
source(file = file.path(path_sourcefile, 'fun_plot_diagnosis.R'))
source(file = file.path(path_sourcefile, 'fun_plot_heatmap.R'))
source(file = file.path(path_sourcefile, 'fun_pipeline_main.R'))

## use parameters 
## PennCNV ( sample LRR mean and SD )
samples_LRR <- read.delim(file = file.path(path_data, "samples_QC.txt"), as.is = TRUE)
samples_LRR$Sample_ID <- sub("\\.txt$", "", samples_LRR$File)

## dup pairs with column_name ( sample1.name sample2.name )
dup_pairs <- read.delim(file = file.path(path_data, "duplicate_pairs.txt"), as.is = TRUE)

## paras_LRR ------------------------------------------------------
paras_LRR <- list(LRR_mean = list(CN_1 = -0.4156184, CN_3 = 0.1734862),
                  LRR_sd   = list(CN_1 = 0.2502591, CN_3 = 0.2249798))  ## sd for one SNP

# main part for runing on cluster -----------------------------------------

cat("read in CNV ...\n")  
dt_cnvs <- read.delim(file = file.path(path_data, "cnv_step3_clean.txt"), as.is = TRUE)

# PennCNV PFB information
cat("read in PFB ...\n")
dt_PFB <- read.table(file = file.path(path_data, "SNP.pfb"), sep = "\t",
                     header = TRUE, as.is = TRUE, check.names = FALSE,
                     comment.char = "")
dt_PFB <- dt_PFB[, c("Name", "PFB", "Position")] # add Position information here

# read in matrix dat of LRR and BAF
cat("read in BAF matrix ...\n")
file_BAF <- paste0("matrix_chr_", chr1, "_BAF.rds")
dt_matrix_BAF <- readRDS(file = file.path( path_matrix, "BAF", file_BAF))
dt_matrix_BAF <- as.matrix(dt_matrix_BAF)

cat("read in LRR matrix ...\n")
file_LRR <- paste0("matrix_chr_", chr1, "_LRR.rds")
dt_matrix_LRR <- readRDS(file = file.path( path_matrix, "LRR", file_LRR))
dt_matrix_LRR <- as.matrix(dt_matrix_LRR)

samples   <- rownames(dt_matrix_LRR)
snps      <- colnames(dt_matrix_LRR)
n_snps    <- length(snps)
n_samples <- length(samples)

# read in cnvrs dat -------------------------------------------------

dt_cnvrs1 <- data.frame()
cnvrs <- NULL  
if (type1 == 0) {
  
  file_cnvr <- "cnvr_batch.txt"  ## with batch information
  dt_cnvrs  <- read.delim(file = file.path(path_data, file_cnvr), as.is = TRUE)
  dt_cnvrs1 <- subset(dt_cnvrs, chr == chr1 & batch == batch1)
  cnvrs <- unique( dt_cnvrs1$CNVR_ID ) 
  
} else if (type1 == 1) {
  
  ## this path can be specified by users
  file_cnvr <- paste0("cnvrs_chr_", chr1, "_batch_", batch1, "_failed.rds")
  
  dt_cnvrs1 <- readRDS(file = file.path(path_data, file_cnvr))
  cnvrs <- unique( dt_cnvrs1$CNVR_ID ) 
}


create_path <- function(path_main, str_subpath) {
  
  path_sub = file.path(path_main, str_subpath) 
  
  if ( !dir.exists(paths = path_sub) ) {
    dir.create(path = path_sub, showWarnings = FALSE, recursive = TRUE)
  }
  
  return( path_sub )
}

# output pathsub_folder: summary/steps/diag/heatmap
path_main <- path_result
path_log     <- create_path(path_main = path_main, str_subpath = "log")
path_png     <- create_path(path_main = path_main, str_subpath = "png")
path_pred    <- create_path(path_main = path_main, str_subpath = "pred")
path_pars    <- create_path(path_main = path_main, str_subpath = "pars")
path_heatmap <- create_path(path_main = path_png, str_subpath = "heatmap")

# add subfolders for each chr and each batch 
folder.name <- paste0("chr_", chr1, "_batch_", batch1)
path_pred <- file.path(path_pred, folder.name)

# test if folder exist
if ( !dir.exists(path_pred) ) {
  dir.create(path = path_pred, showWarnings = FALSE, recursive = TRUE)
}

## must be changed here to save each CNVRID data
path_cnvr_stat <- file.path(path_result, "stats")
dir.create(path = path_cnvr_stat, showWarnings = FALSE)


res_pars_all <- data.frame() 
# --------------------------------------------------------------------
for (i in 1:nrow(dt_cnvrs1)) {

  cnvr1 <- dt_cnvrs1$CNVR_ID[i]

  cat("cnvr1:", cnvr1, i, "in", nrow(dt_cnvrs1), "\n")

  snp_start <- dt_cnvrs1$start_snp[i]
  snp_end   <- dt_cnvrs1$end_snp[i]

  ## snps have been sorted by their positions on chromosome
  ## when preparing chromosome-wise LRR and BAF matrices 
  idx_start <- which(snps == snp_start)
  idx_end   <- which(snps == snp_end)
  
  ## check idx_start and idx_end
  idxs <- c(idx_start, idx_end)
  if (length(idxs) != 2 | idx_start >= idx_end) {
    stop("CNVR boundaries are not consistency with SNP information.")
  }

  snps_name <- snps[idx_start:idx_end] # all snps in cnvr1

  # plot heatmap add 20 snps on the both side ----------------------------
  idx_outer_start <- dt_cnvrs1$outer.start[i]
  idx_outer_end   <- dt_cnvrs1$outer.end[i]
  idx_start_new <- ifelse((idx_outer_start - 20) <= 0, 1, idx_outer_start - 20)  # new start and end for plot heatmap
  idx_end_new <- ifelse((idx_outer_end + 20) > n_snps, n_snps, idx_outer_end + 20)

  dt_lrr_heatmap = dt_matrix_LRR[, idx_start_new:idx_end_new]

  snps_name_heatmap <- snps[idx_start_new:idx_end_new]
  snps_name_all     <- snps[idx_outer_start:idx_outer_end]
  # colnames(dt_lrr_heatmap) <- snps_name_heatmap

  snps_add   <- setdiff(snps_name_heatmap, snps_name_all)
  snps_outer <- setdiff(snps_name_all, snps_name)
  snps_flag  <- ifelse(snps_name_heatmap %in% snps_add, 0,
                      ifelse(snps_name_heatmap %in% snps_name, 2, 1))
  dt_snps_flag <- data.frame(snp_name = snps_name_heatmap,
                             snp_flag = snps_flag,
                             stringsAsFactors = FALSE)

  filename_heatmap <- paste0("heatmap_", cnvr1, ".png")
  png(filename = file.path(file.path(path_png, "heatmap"), filename_heatmap),
      width = 12, height = 12, units = "in", res = 512)
  plot_heatmap(dt_lrr_heatmap = dt_lrr_heatmap, dt_snps_flag = dt_snps_flag)
  dev.off()
  
  # -------------------------------------------------------------------
  dt_baf = dt_matrix_BAF[, idx_start:idx_end]
  dt_lrr = dt_matrix_LRR[, idx_start:idx_end]
  
  numsnp <- idx_end - idx_start + 1
  samples_new <- rownames(dt_baf) ## need change in dt_cnvr_stat
  stopifnot( all(samples_new == samples) )

  dt_cnvr_stat <- data.frame(CNVR_ID = cnvr1, 
                             Chr = chr1, 
                             BAF = as.vector(dt_baf),
                             LRR = as.vector(dt_lrr), 
                             Sample_ID = rep(samples_new, numsnp),
                             Name = rep(snps_name, each = length(samples_new)), 
                             numSNP = numsnp,
                             stringsAsFactors = FALSE)

  dt_PFB1 <- subset(dt_PFB, Name %in% snps_name)
  dt_cnvr_stat <- merge(dt_cnvr_stat, dt_PFB1, all.x = TRUE)

  dt_samples_cn <- data.frame(Sample_ID = samples_new, stringsAsFactors = FALSE)

  dt_cnv <- subset(dt_cnvs, CNVR_ID == cnvr1)
  dt_cnv <- dt_cnv[, c("Sample_ID", "CN", "alg")]
  dt_samples_cn <- merge(dt_samples_cn, dt_cnv, all.x = TRUE)
  dt_samples_cn$CN[ which(is.na(dt_samples_cn$CN)) ]   <- 2
  dt_samples_cn$alg[ which(is.na(dt_samples_cn$alg)) ] <- "other"

  dt_cnvr_stat <- merge(dt_cnvr_stat, dt_samples_cn, all.x = TRUE)
  
  ## save CNVR-stat data
  saveRDS(dt_cnvr_stat, file = file.path(path_cnvr_stat, paste0(cnvr1, "_stat.rds")))
  
  res_pipeline_cnvr1 <- pipeline_main(dt_cnvrs = dt_cnvr_stat, 
                                      paras_LRR = paras_LRR, 
                                      dup_pairs = dup_pairs,
                                      samples_LRR = samples_LRR, 
                                      path_png = path_png,
                                      n.sample = n_samples)

  res_gatk_pred_final <- res_pipeline_cnvr1$res_gatk_pred_final
  res_pars <- res_pipeline_cnvr1$res_pars

  cat( names(res_pars_all), "\n")
  cat( names(res_gatk_pred_final), "\n")
  
  # res_pred_all <- rbind(res_pred_all, res_gatk_pred_final)
  filename_cnvr1 <- paste0(cnvr1, "_pred.rds")
  saveRDS(res_gatk_pred_final, file = file.path(path_pred, filename_cnvr1))

  res_pars_all <- rbind(res_pars_all, res_pars)

}

filename_pars <- paste0("CNVR_pars_chr_", chr1, "_batch_", batch1, ".rds")
saveRDS(res_pars_all, file = file.path(path_pars, filename_pars)) # pars file





