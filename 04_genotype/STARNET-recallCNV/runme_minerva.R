#!/hpc/packages/minerva-common/R/3.3.1/lib64/R/bin/Rscript --vanilla

args <- commandArgs(trailingOnly = TRUE)

chr1 <- as.integer(args[1])   # chr infor
batch1 <- as.integer(args[2]) # which batch

# load packages
suppressMessages(library(dplyr))
suppressMessages(library(mixtools))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(plyr))
suppressMessages(library(modeest))
suppressMessages(library(mclust))
suppressMessages(library(gridExtra))
suppressMessages(library(pheatmap))
suppressMessages(library(RColorBrewer))

path_sourcefile <- "/sc/orga/projects/haok01a/chengh04/StarNet/work_3/Predict_All_CNVR/model_test_numSNP/functions_gatk"
source(file = file.path(path_sourcefile, 'fun_BAF.R'))
source(file = file.path(path_sourcefile, 'fun_gatk.R'))
source(file = file.path(path_sourcefile, 'fun_LRR.R'))
source(file = file.path(path_sourcefile, 'fun_models_zz_8_cc_raw.R'))
source(file = file.path(path_sourcefile, 'fun_plot_steps.R'))
source(file = file.path(path_sourcefile, 'fun_plot_diagnosis.R'))
source(file = file.path(path_sourcefile, 'fun_plot_heatmap.R'))
source(file = file.path(path_sourcefile, 'fun_pipeline_main.R'))

# use parameters 
# parameters
## PennCNV sample LRR
samples_LRR <- read.table('/sc/orga/projects/haok01a/chengh04/StarNet/work_3/Predict_All_CNVR/res_one_CNVR/StarNet_PennCNV_qc.txt', sep = "\t",
                          header = TRUE, as.is = TRUE, check.names = FALSE)
samples_LRR$Sample_ID <- gsub(pattern = ".txt", replacement = "", samples_LRR$File)

## sample1.name sample2.name
dup_pairs <- read.table('/sc/orga/projects/haok01a/chengh04/StarNet/work_3/Predict_All_CNVR/res_one_CNVR/samples_dup_pair_qc.txt', sep = "\t",
                        header = TRUE, as.is = TRUE, check.names = FALSE)

## paras_LRR
paras_LRR <- list(LRR_mean = list(CN_1 = -0.4156184, CN_3 = 0.1734862),
                  LRR_sd   = list(CN_1 = 0.2502591, CN_3 = 0.2249798))  # sd for one SNP

# --- main for runing on minerva -----
file_cnvr <- "cnvrs_boundery_with_batch_info.rds"

path_cnvr <- "/sc/orga/projects/haok01a/chengh04/StarNet/work_3/Predict_All_CNVR/model_test_numSNP/data"
# all samples
samples <- readRDS(file = file.path(path_cnvr, "samples.rds"))  # all 834 samples in STARNET
# all cnvs 
dt_cnvs <- readRDS(file = file.path(path_cnvr, "cnvs_UnifiedCNV_new.rds"))
# PFB 
dt_PFB <- read.delim(file = file.path(path_cnvr, "SNP.pfb"), as.is = TRUE)
dt_PFB <- dt_PFB[, c("Name", "PFB", "Position")] # add Position information here

# read in matrix dat of LRR and BAF
path_matrix <- "/sc/orga/projects/haok01a/chengh04/StarNet/data/matrix_chr_based_LRR_BAF"
file_BAF <- paste0("matrix_chr_", chr1, "_BAF.rds")
dt_matrix_BAF <- readRDS(file = file.path(paste(path_matrix, "BAF", sep = "/"), file_BAF))

samples_all <- colnames(dt_matrix_BAF)
idxs_sub_col <- which(samples_all %in% samples)

dt_matrix_BAF <- dt_matrix_BAF[, idxs_sub_col]  # matrix BAF used

file_LRR <- paste0("matrix_chr_", chr1, "_LRR.rds")
dt_matrix_LRR <- readRDS(file = file.path(paste(path_matrix, "LRR", sep = "/"), file_LRR))

dt_matrix_LRR <- dt_matrix_LRR[, idxs_sub_col]  # matrix LRR used

file_annotate <- paste0("chr_", chr1, "_annotate.rds")  # annotation information
dt_annotate <- readRDS(file = file.path(paste(path_matrix, "annotate", sep = "/"), file_annotate))

# read in cnvrs dat -------------------------------------------------
dt_cnvrs  <- readRDS(file = file.path(path_cnvr, file_cnvr))
dt_cnvrs1 <- subset(dt_cnvrs, chr == chr1 & batch == batch1)
cnvrs <- unique(dt_cnvrs1$CNVR_ID)  ## all cnvrs

# output path
# sub folder: summary/steps/diag/heatmap
path_png<- "/sc/orga/projects/haok01a/chengh04/StarNet/work_3/Predict_All_CNVR/model_test_numSNP/functions_gatk/res/png"
path_pred <- "/sc/orga/projects/haok01a/chengh04/StarNet/work_3/Predict_All_CNVR/model_test_numSNP/functions_gatk/res/pred"
path_pars <- "/sc/orga/projects/haok01a/chengh04/StarNet/work_3/Predict_All_CNVR/model_test_numSNP/functions_gatk/res/pars"

# res_pred_all <- data.frame() ## save each CNVR
res_pars_all <- data.frame() ##
for (i in 1:nrow(dt_cnvrs1)) {
  
  cnvr_id <- dt_cnvrs1$CNVR_ID[i]
  
  cat("CNVR_ID:", cnvr_id, i, "in", nrow(dt_cnvrs1), "\n")
  
  snp_start <- dt_cnvrs1$start_snp[i]
  snp_end   <- dt_cnvrs1$end_snp[i]
  
  idx_start <- which(dt_annotate$SNP_Name == snp_start)
  idx_end   <- which(dt_annotate$SNP_Name == snp_end)
  
  snps_name <- dt_annotate$SNP_Name[idx_start:idx_end] # all snps name
  
  # plot heatmap
  # add 20 snps on the both side --------------------------------------------
  idx_outer_start <- dt_cnvrs1$outer.start[i]
  idx_outer_end   <- dt_cnvrs1$outer.end[i] 
  idx_start_new <- ifelse(idx_outer_start - 20 <= 0, 1, idx_outer_start - 20)  # new start and end for plot heatmap
  idx_end_new <- ifelse(idx_outer_end + 20 > nrow(dt_annotate), nrow(dt_annotate), idx_outer_end + 20)
  
  dt_lrr_heatmap <- dt_matrix_LRR[idx_start_new:idx_end_new, ]
  dt_lrr_heatmap <- t(dt_lrr_heatmap)  ## transform
  
  snps_name_heatmap <- dt_annotate$SNP_Name[idx_start_new:idx_end_new]
  snps_name_all <- dt_annotate$SNP_Name[idx_outer_start:idx_outer_end]
  colnames(dt_lrr_heatmap) <- snps_name_heatmap
  
  snps_add <- setdiff(snps_name_heatmap, snps_name_all)
  snps_outer <- setdiff(snps_name_all, snps_name)
  snps_flag <- ifelse(snps_name_heatmap %in% snps_add, 0, 
                      ifelse(snps_name_heatmap %in% snps_name, 2, 1))
  dt_snps_flag <- data.frame(snp_name = snps_name_heatmap, 
                             snp_flag = snps_flag, 
                             stringsAsFactors = FALSE)
  
  filename_heatmap <- paste0("heatmap_", cnvr_id, ".png")
  png(filename = file.path(file.path(path_png, "heatmap"), filename_heatmap), width = 12, height = 12, units = "in", res = 512)
  plot_heatmap(dt_lrr_heatmap = dt_lrr_heatmap, dt_snps_flag = dt_snps_flag)
  dev.off()
  # --------------------------------------------------------------------
  
  idxs <- c(idx_start, idx_end)
  if (length(idxs) != 2 | idx_start >= idx_end) {
    stop("Not consistency with SNP information.")
  }
  
  dt_baf <- dt_matrix_BAF[idx_start:idx_end, ]
  dt_lrr <- dt_matrix_LRR[idx_start:idx_end, ]  
  
  numsnp <- idx_end - idx_start + 1
  samples_new <- colnames(dt_baf) ## need change in dt_cnvr_stat
  
  dt_cnvr_stat <- data.frame(CNVR_ID = cnvr_id, Chr = chr1, BAF = as.vector(dt_baf), 
                             LRR = as.vector(dt_lrr), Sample_ID = rep(samples_new, each = numsnp),
                             Name = rep(snps_name, length(samples_new)), numSNP = numsnp, 
                             stringsAsFactors = FALSE)
  
  dt_PFB1 <- subset(dt_PFB, Name %in% snps_name)
  dt_cnvr_stat <- merge(dt_cnvr_stat, dt_PFB1, all.x = TRUE)
  
  dt_samples_cn <- data.frame(Sample_ID = samples, stringsAsFactors = FALSE)
  
  dt_cnv <- subset(dt_cnvs, CNVR_ID == cnvr_id) #
  dt_cnv <- dt_cnv[, c("Sample_ID", "CN", "alg")]
  dt_samples_cn <- merge(dt_samples_cn, dt_cnv, all.x = TRUE)
  dt_samples_cn$CN[which(is.na(dt_samples_cn$CN))]   <- 2
  dt_samples_cn$alg[which(is.na(dt_samples_cn$alg))] <- "other"
  
  dt_cnvr_stat <- merge(dt_cnvr_stat, dt_samples_cn, all.x = TRUE)  ##
  
  res_pipeline_cnvr1 <- pipeline_main(dt_cnvrs = dt_cnvr_stat, paras_LRR = paras_LRR, dup_pairs = dup_pairs,
                                      samples_LRR = samples_LRR, path_png = path_png)
  
  res_gatk_pred_final <- res_pipeline_cnvr1$res_gatk_pred_final
  res_pars <- res_pipeline_cnvr1$res_pars
  
  # res_pred_all <- rbind(res_pred_all, res_gatk_pred_final)
  filename_cnvr1 <- paste0(cnvr_id, "_pred.rds")
  saveRDS(res_gatk_pred_final, 
          file = file.path(path_pred, filename_cnvr1))
  
  res_pars_all <- rbind(res_pars_all, res_pars)
  
}

# filename_pred <- paste0("CNVR_pred_chr_", chr1, "_batch_", batch1, ".rds")
filename_pars <- paste0("CNVR_pars_chr_", chr1, "_batch_", batch1, ".rds")

# saveRDS(res_pred_all, file = file.path(path_pred, filename_pred)) # pred file
saveRDS(res_pars_all, file = file.path(path_pars, filename_pars)) # pars file 




