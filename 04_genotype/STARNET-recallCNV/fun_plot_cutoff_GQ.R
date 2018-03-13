#!/hpc/packages/minerva-common/R/3.3.1/lib64/R/bin/Rscript --vanilla

suppressMessages(library(plyr))
## sample1.name sample2.name --------------------------------------------------
dup_pairs <- read.table('/sc/orga/projects/haok01a/chengh04/StarNet/work_3/Predict_All_CNVR/res_one_CNVR/samples_dup_pair_qc.txt', sep = "\t",
                        header = TRUE, as.is = TRUE, check.names = FALSE)
n_dup_pair <- nrow(dup_pairs) # number of duplicate pairs

file_cnvr <- "cnvrs_boundery_with_batch_info.rds"
path_cnvr <- "/sc/orga/projects/haok01a/chengh04/StarNet/work_3/Predict_All_CNVR/model_test_numSNP/data"
dt_cnvrs  <- readRDS(file = file.path(path_cnvr, file_cnvr))
n_cnvr    <- nrow(dt_cnvrs)
n_sample  <- 834  # need to be changed

path_cnvr_pred <- "/sc/orga/projects/haok01a/chengh04/StarNet/work_3/Predict_All_CNVR/model_test_numSNP/functions_gatk/res/pred"

# filenames <- list.files(path = path_cnvr_pred)  # list all cnvr_prediction
# # combine ------------------------------------------------------------
# ns <- n_cnvr*n_sample
# res_pred_simple <- data.frame(Sample_ID = character(ns),
#                               CNVR_ID = character(ns),
#                               value_GQ = numeric(ns),
#                               CN_gatk_pred = integer(ns), stringsAsFactors = FALSE)
# for (i in 1:length(filenames)) {
#   
#   f1 <- filenames[i]
#   cat("Read in:", f1, i, "in", length(filenames), "\n")
#   dat1 <- readRDS(file = file.path(path_cnvr_pred, f1))
#   dat1 <- dat1[, c("Sample_ID", "CNVR_ID", "value_GQ", "CN_gatk_pred")]
#   
#   idxs <- ((i-1)*n_sample+1):(n_sample*i)
#   res_pred_simple[idxs, ] <- dat1 
#   
# }
# 
# path_output <- "/sc/orga/projects/haok01a/chengh04/StarNet/work_3/Predict_All_CNVR/model_test_numSNP/functions_gatk/res/res_summary"
# saveRDS(res_pred_simple, file = file.path(path_output, "dat_cnvr_all_pred.rds"))

# read in data ---------------------------------------------------
path_output <- "/sc/orga/projects/haok01a/chengh04/StarNet/work_3/Predict_All_CNVR/model_test_numSNP/functions_gatk/res/res_summary"
res_pred_simple <- readRDS(file = file.path(path_output, "dat_cnvr_all_pred.rds"))

n_sample <- 834
n_cnvr <- 4006

# function ------------------------------------------------
calculate_dup_pair_consistency <- function(dat, dup_pairs) {
  
  # clean cnvr with all CN = 2 and CN = 4 
  dat <- dat[order(dat$Sample, dat$CNVR_ID), ]
  samples <- unique(dat$Sample)
  cnvrs   <- dat$CNVR_ID[1:4006]
  mat <- matrix(dat$CN_gatk_pred, nrow = 4006, ncol = 834,
                dimnames = list(cnvrs, samples))
  
  # filter CNVR
  freqs <- unlist(lapply(1:4006, FUN = function(k) {
    v1 <- as.vector(mat[k, ])
    freq1 <- sum(v1 %in% c(0, 1, 3))
    freq1
  }))
  
  idxs_del <- which(freqs == 0)
  mat1 <- mat[-idxs_del, ]
  n_cnvr <- nrow(mat1)
  
  consistency_rates <- c()
  for (i in 1:nrow(dup_pairs)) {
    
    sample1 <- dup_pairs$sample1.name[i]
    sample2 <- dup_pairs$sample2.name[i]
    
    cns1 <- as.vector(mat1[, sample1])
    cns2 <- as.vector(mat1[, sample2])   # copy number of sample2
    
    # filter nocall cnvr
    idxs <- union(which(cns1 == 4), which(cns2 == 4))
    if(length(idxs) >= 1) {
      cns1 <- cns1[-idxs]
      cns2 <- cns2[-idxs]
    }
    
    idxs_overlap <- which(cns1 != 2 & cns2 != 2 & cns1 == cns2)
    idxs_union   <- union(which(cns1 != 2), which(cns2 != 2))
    
    res1 <- length(idxs_overlap)/length(idxs_union)
    consistency_rates  <- c(consistency_rates, res1)
  }
  
  res <- data.frame(consistency_rate = consistency_rates,
                    sample1.name = dup_pairs$sample1.name,
                    sample2.name = dup_pairs$sample2.name,
                    n_cnvr = n_cnvr,
                    stringsAsFactors = FALSE)

  res  ## return
}

## clean matrix
generate_matrix_clean <- function(dat) {
  
  dat <- dat[order(dat$Sample_ID, dat$CNVR_ID), ]
  samples <- unique(dat$Sample_ID)
  cnvrs <- dat$CNVR_ID[1:4006]
  
  mat <- matrix(dat$CN_gatk_pred, nrow = length(cnvrs),
                ncol = length(samples),
                dimnames = list(cnvrs, samples))
  
  freqs <- unlist(lapply(1:nrow(mat), FUN = function(k) {
    v1 <- as.vector(mat[k, ])
    freq1 <- sum(v1 %in% c(0, 1, 3))
    freq1
  }))
  
  idxs_del <- which(freqs == 0)
  if (length(idxs_del) >= 0) {
    mat1 <- mat[-idxs_del, ]
    return(mat1)
  } else {
    return(mat)
  }
  
}

# cutoff of GQ -------------------------------------------------------
cutoffs_gq <- c(0, 5, 10, 15, 20, 30, 50, 60, 80, 100) ## change
# call_rate and boxplot(dup pair consistency)
# value_GQ/CN_gatk_pred
res_dup_pair_consistency <- data.frame()
res_call_rate_sample <- data.frame()
res_call_rate_cnvr <- data.frame()
res_barplot <- data.frame()
for (i in 1:length(cutoffs_gq)) {
  
  gq1 <- cutoffs_gq[i]
  idxs_nocall <- which(res_pred_simple$value_GQ <= gq1)
  
  res_pred_simple_gq <- res_pred_simple
  if (length(idxs_nocall) >= 1) {
    res_pred_simple_gq$CN_gatk_pred[idxs_nocall] <- 4 ## nocall set as 4
  }
  
  # consistency rate
  dups_consistency <- calculate_dup_pair_consistency(dat = res_pred_simple_gq, dup_pairs = dup_pairs)
  # res_dups <- data.frame(cutoff_gq = rep(gq1, n_dup_pair), 
  #                        consistency_rate = dups_consistency)
  dups_consistency$cutoff_gq <- gq1
  cat("median:", median(dups_consistency$consistency_rate), "\n")
  
  res_dup_pair_consistency <- rbind(res_dup_pair_consistency, dups_consistency)  ## for boxplot
  
  # res_sample <- res_pred_simple_gq[order(res_pred_simple_gq$Sample_ID), ]
  # samples <- unique(res_sample$Sample_ID)
  mat_clean <- generate_matrix_clean(dat = res_pred_simple_gq)
  n_sample_clean <- ncol(mat_clean)
  n_cnvr_clean <- nrow(mat_clean)
  
  cat("cutoff of GQ score:", gq1, "CNVR:", n_cnvr_clean, "Samples:", n_sample_clean, "\n") # print out
  call_rate_sample <- unlist(lapply(1:n_sample_clean, FUN = function(k) {
    v1 <- mat_clean[, k]
    f1 <- sum(v1 %in% c(0, 1, 2, 3))
    f1/n_cnvr_clean
  }))
  res_sample_gq1 <- data.frame(cutoff_gq = rep(gq1, n_sample_clean), 
                               Sample_ID = colnames(mat_clean),
                               call_rate = call_rate_sample, stringsAsFactors = FALSE)
  
  call_rate_cnvr <- unlist(lapply(1:n_cnvr_clean, FUN = function(k) {
    v1 <- mat_clean[k, ]
    f1 <- sum(v1 %in% c(0, 1, 2, 3))
    f1/n_sample_clean  # call rate
  }))
  res_cnvr_gq1 <- data.frame(cutoff_gq = rep(gq1, n_cnvr_clean),
                             CNVR_ID = rownames(mat_clean),
                             call_rate = call_rate_cnvr, stringsAsFactors = FALSE)
  
  res_call_rate_cnvr <- rbind(res_call_rate_cnvr, res_cnvr_gq1)
  res_call_rate_sample <- rbind(res_call_rate_sample, res_sample_gq1)
  
  res_bar1 <- data.frame(GQ_score = gq1, number_CNVR = n_cnvr_clean)
  res_barplot <- rbind(res_barplot, res_bar1)
}

filename_dup_pair_consistency <- "dat_dup_pair_consistency.rds"
filename_call_rate_sample <- "dat_call_rate_sample.rds"
filename_call_rate_cnvr <- "dat_call_rate_cnvr.rds"

saveRDS(res_dup_pair_consistency, file = file.path(path_output, filename_dup_pair_consistency))
saveRDS(res_call_rate_cnvr, file = file.path(path_output, filename_call_rate_cnvr))
saveRDS(res_call_rate_sample, file = file.path(path_output, filename_call_rate_sample))


# plot
res_dup_pair_consistency$cutoff_gq <- factor(res_dup_pair_consistency$cutoff_gq, 
                                             levels = cutoffs_gq)
png(filename = file.path(path_output, "eCNV_dup_pair_consistent_rate_as_a_function_of_GQ_score_cutoff.png"),
    width = 12, height = 12, units = "in", res = 512)
boxplot(consistency_rate ~ cutoff_gq, data = res_dup_pair_consistency,
        xlab = "cutoff of GQ score", ylab = "consistency rate")
dev.off()


res_call_rate_cnvr$cutoff_gq <- factor(res_call_rate_cnvr$cutoff_gq, levels = cutoffs_gq)
png(filename = file.path(path_output, "boxplot_medians_cnvr_based.png"),
    width = 12, height = 12, units = "in", res = 512)
boxplot(call_rate ~ cutoff_gq, data = res_call_rate_cnvr, 
        xlab = "cutoff of GQ score", ylab = "call_rate", 
        main = "median of call_rate cnvr based")
dev.off()

res_call_rate_sample$cutoff_gq <- factor(res_call_rate_sample$cutoff_gq, levels = cutoffs_gq)
png(filename = file.path(path_output, "boxplot_medians_sample_based.png"),
    width = 12, height = 12, units = "in", res = 512)
boxplot(call_rate ~ cutoff_gq, data = res_call_rate_sample, 
        xlab = "cutoff of GQ score", ylab = "call_rate", 
        main = "median of call_rate sample based")
dev.off()


# barplot 
# number of CNVR in each QG cutoff
png(filename = file.path(path_output, "number_cnvr_in_each_GQ.png"),
    width = 12, height = 12, units = "in", res = 512)

res_barplot$GQ_score <- factor(res_barplot$GQ_score, levels = cutoffs_gq)
ggplot(data = res_barplot, aes(GQ_score, number_CNVR)) + 
  geom_col() + 
  theme_bw(base_size = 9) + 
  labs(x = "GQ score", y = "number of CNVR")

dev.off()






# res_call_rate_cnvr$cutoff_gq <- factor(res_call_rate_cnvr$cutoff_gq, 
#                                        levels = cutoffs_gq)
# dat1 <- ddply(res_call_rate_cnvr, .(cutoff_gq), summarize, median = median(call_rate))
# dat1$cutoff_gq <- as.integer(as.character(dat1$cutoff_gq))
# png(filename = file.path(path_output, "medians_cnvr_based.png"),
#     width = 12, height = 12, units = "in", res = 512)
# plot(median~cutoff_gq, data = dat1, type = "b",
#      xlab = "cutoff of QG score", ylab = "median of nocall_rate",
#      main = "median of nocall_rate cnvr based")
# dev.off()
# 
# res_call_rate_sample$cutoff_gq <- factor(res_call_rate_sample$cutoff_gq,
#                                          levels = cutoffs_gq)
# dat2 <- ddply(res_call_rate_sample, .(cutoff_gq), summarize, median = median(call_rate))
# dat2$cutoff_gq <- as.integer(as.character(dat1$cutoff_gq))
# png(filename = file.path(path_output, "medians_sample_based.png"),
#     width = 12, height = 12, units = "in", res = 512)
# plot(median~cutoff_gq, data = dat2, type = "b",
#      xlab = "cutoff of QG score", ylab = "median of nocall_rate",
#      main = "median of nocall_rate sample based")
# dev.off()





