
library(plyr)
path_output <- "/sc/orga/projects/haok01a/chengh04/StarNet/work_3/Predict_All_CNVR/model_test_numSNP/functions_gatk/res/res_summary"
res_pred_simple <- readRDS(file = file.path(path_output, "dat_cnvr_all_pred.rds"))
n_sample = 834
n_cnvr = 4006

# cutoff of GQ -------------------------------------------------------
cutoffs_gq <- c(0, 5, 10, 15, 20, 30, 40, 50, 80, 100) ## change
# call_rate and boxplot(dup pair consistency)
# value_GQ/CN_gatk_pred
res_call_rate_sample <- data.frame()
res_call_rate_cnvr <- data.frame()

for (i in 1:length(cutoffs_gq)) {
  
  gq1 <- cutoffs_gq[i]
  cat("cutoff of GQ score:", gq1, "\n") # print out
  idxs_nocall <- which(res_pred_simple$value_GQ <= gq1)
  
  res_pred_simple_gq <- res_pred_simple
  if (length(idxs_nocall) >= 1) {
    res_pred_simple_gq$CN_gatk_pred[idxs_nocall] <- 4 ## nocall set as 4
  }
  
  res_sample <- res_pred_simple_gq[order(res_pred_simple_gq$Sample_ID), ]
  samples <- unique(res_sample$Sample_ID)
  call_rate_sample <- unlist(lapply(1:n_sample, FUN = function(k) {
    idxs <- ((k-1)*n_cnvr+1):(n_cnvr*k)
    cns <- res_sample$CN_gatk_pred[idxs]
    length(which(cns != 4))/length(cns) # call_rate
  }))
  res_sample_gq1 <- data.frame(cutoff_gq = rep(gq1, n_sample), 
                               Sample_ID = samples,
                               call_rate = call_rate_sample, stringsAsFactors = FALSE)
  
  res_cnvr <- res_pred_simple_gq[order(res_pred_simple_gq$CNVR_ID), ]
  cnvrs <- unique(res_cnvr$CNVR_ID)
  call_rate_cnvr <- unlist(lapply(1:n_cnvr, FUN = function(k) {
    idxs <- ((k-1)*n_sample+1):(n_sample*k)
    cns <- res_cnvr$CN_gatk_pred[idxs]
    length(which(cns != 4))/length(cns)  # call rate
  }))
  res_cnvr_gq1 <- data.frame(cutoff_gq = rep(gq1, n_cnvr),
                             CNVR_ID = cnvrs,
                             call_rate = call_rate_cnvr, stringsAsFactors = FALSE)
  
  res_call_rate_cnvr <- rbind(res_call_rate_cnvr, res_cnvr_gq1)
  res_call_rate_sample <- rbind(res_call_rate_sample, res_sample_gq1)
}


# boxplot
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

# CNVR_67_r1_chr3_q
# CNVR_54_r1_chr7_q
# CNVR_68_r1_chr6_q

# median scatter plot
res_call_rate_cnvr$cutoff_gq <- factor(res_call_rate_cnvr$cutoff_gq, 
                                       levels = cutoffs_gq)
dat1 <- ddply(res_call_rate_cnvr, .(cutoff_gq), summarize, median = median(call_rate))
dat1$cutoff_gq <- as.integer(as.character(dat1$cutoff_gq))
png(filename = file.path(path_output, "medians_cnvr_based.png"),
    width = 12, height = 12, units = "in", res = 512)
plot(median~cutoff_gq, data = dat1, type = "b",
     xlab = "cutoff of QG score", ylab = "median of call_rate",
     main = "median of call_rate cnvr based")
dev.off()

res_call_rate_sample$cutoff_gq <- factor(res_call_rate_sample$cutoff_gq,
                                         levels = cutoffs_gq)
dat2 <- ddply(res_call_rate_sample, .(cutoff_gq), summarize, median = median(call_rate))
dat2$cutoff_gq <- as.integer(as.character(dat1$cutoff_gq))
png(filename = file.path(path_output, "medians_sample_based.png"),
    width = 12, height = 12, units = "in", res = 512)
plot(median~cutoff_gq, data = dat2, type = "b",
     xlab = "cutoff of QG score", ylab = "median of call_rate",
     main = "median of call_rate sample based")
dev.off()

