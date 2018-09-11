
cnvr1 <- "cnvr1"
chr1 <- 1
path_to_script <- "./script"
path_to_data <- "./test_data_regenotype"
path_to_res  <- "./res"

## this script for regenotype CNVR (only on CNVR)

suppressMessages({
  library(dplyr)
  library(mixtools)
  library(ggplot2)
  library(cowplot)
  library(plyr)
  library(modeest)
  library(mclust)
  library(gridExtra)
  library(pheatmap)
  library(RColorBrewer)
})

## load functions in script folder
path_sourcefile <- path_to_script
source(file = file.path(path_sourcefile, 'fun_BAF.R'))
source(file = file.path(path_sourcefile, 'fun_gatk.R'))
source(file = file.path(path_sourcefile, 'fun_LRR.R'))
source(file = file.path(path_sourcefile, 'fun_models_STARNET.R'))
source(file = file.path(path_sourcefile, 'fun_plot_steps.R'))
source(file = file.path(path_sourcefile, 'fun_plot_diagnosis.R'))
source(file = file.path(path_sourcefile, 'fun_plot_heatmap.R'))
source(file = file.path(path_sourcefile, 'fun_pipeline_main.R'))

## paras_LRR
paras_LRR <- list(LRR_mean = list(CN_1 = -0.4156184, CN_3 = 0.1734862),
                  LRR_sd   = list(CN_1 = 0.2502591, CN_3 = 0.2249798))  # sd for one SNP

## put all input data in data folder
path_dat <- path_to_data

samples_LRR <- read.table(file = file.path(path_dat, 'StarNet_PennCNV_qc.txt'), sep = "\t",
                          header = TRUE, as.is = TRUE, check.names = FALSE)
samples_LRR$Sample_ID <- gsub(pattern = ".txt", replacement = "", samples_LRR$File)

dup_pairs <- read.table(file = file.path(path_dat, 'samples_dup_pair_qc.txt'), sep = "\t",
                        header = TRUE, as.is = TRUE, check.names = FALSE) ## sample1.name sample2.name

dat_cnvr1 <- readRDS(file = file.path(path_dat, 'dat_cnvr1.rds'))

dt_matrix_BAF <- readRDS(file = file.path(path_dat, "matrix_BAF_cnvr1.rds"))
dt_matrix_LRR <- readRDS(file = file.path(path_dat, "matrix_LRR_cnvr1.rds"))

samples <- rownames(dt_matrix_BAF)
snps    <- colnames(dt_matrix_BAF)

dat_PFB_chr1 <- read.delim(file = file.path(path_dat, "SNP.cnvr1.pfb"), as.is = TRUE)
dat_PFB_chr1 <- dat_PFB_chr1[, c("Name", "Chr", "PFB", "Position")] # add Position information here
nrow( dat_PFB_chr1 )

## save all output results in res folder
path_res <- path_to_res
dat_cnvs1 <- readRDS(file = file.path(path_res, "cnvs_step3_clean.rds")) # from test_ensembleCNV step

# create subfolder
create_path <- function(path_main, str_subpath) {
  path_sub = file.path(path_main, str_subpath) 
  if ( !dir.exists(paths = path_sub) ) {
    dir.create(path = path_sub, showWarnings = FALSE, recursive = TRUE)
  }
  return( path_sub )
}

path_log  = create_path(path_main = path_res, str_subpath = "log")
path_png  = create_path(path_main = path_res, str_subpath = "png")
path_pred = create_path(path_main = path_res, str_subpath = "pred")
path_pars = create_path(path_main = path_res, str_subpath = "pars")

# path_heatmap = create_path(path_main = path_png, str_subpath = "heatmap")

## start training
snp_start <- dat_cnvr1$start_snp
snp_end   <- dat_cnvr1$end_snp

idx_start <- which(snps == snp_start)
idx_end   <- which(snps == snp_end)
  
## check idx_start and idx_end
idxs <- c(idx_start, idx_end)
if (length(idxs) != 2 | idx_start >= idx_end) {
  stop("Not consistency with snpName information.")
}

snps_name <- snps[idx_start:idx_end] # all snps in cnvr1

dt_baf = dt_matrix_BAF[, idx_start:idx_end]
dt_lrr = dt_matrix_LRR[, idx_start:idx_end]
  
numsnp <- idx_end - idx_start + 1
samples_new <- rownames(dt_baf) ## need change in dt_cnvr_stat
stopifnot( all(samples_new == samples) )

n.sample <- length( samples_new ) 

dt_cnvr_stat <- data.frame(CNVR_ID = cnvr1, Chr = chr1, BAF = as.vector(dt_baf),
                           LRR = as.vector(dt_lrr), Sample_ID = rep(samples_new, numsnp),
                           Name = rep(snps_name, each = length(samples_new)), numSNP = numsnp,
                           stringsAsFactors = FALSE)

dt_cnvr_stat <- merge(dt_cnvr_stat, dt_PFB1, all.x = TRUE)
dt_samples_cn <- data.frame(Sample_ID = samples_new, stringsAsFactors = FALSE)

dat_cnvs1 <- dat_cnvs1[, c("Sample_ID", "CN", "alg")]
dt_samples_cn <- merge(dt_samples_cn, dat_cnvs1, all.x = TRUE)
dt_samples_cn$CN[which(is.na(dt_samples_cn$CN))]   <- 2
dt_samples_cn$alg[which(is.na(dt_samples_cn$alg))] <- "other"

dt_cnvr_stat <- merge(dt_cnvr_stat, dt_samples_cn, all.x = TRUE)  ##

res_pipeline_cnvr1 <- pipeline_main(dt_cnvrs = dt_cnvr_stat, paras_LRR = paras_LRR, dup_pairs = dup_pairs,
                                    samples_LRR = samples_LRR, path_png = path_png, n.sample = n.sample)

res_gatk_pred_final <- res_pipeline_cnvr1$res_gatk_pred_final
res_pars <- res_pipeline_cnvr1$res_pars

saveRDS(res_gatk_pred_final, file = file.path(path_pred, paste0(cnvr1, ".pred.rds")))
saveRDS(res_pars, file = file.path(path_pars, paste0(cnvr1, ".pars.rds")))






