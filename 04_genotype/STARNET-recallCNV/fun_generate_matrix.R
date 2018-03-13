
# output path named: res_summary
path_output <- "/sc/orga/projects/haok01a/chengh04/StarNet/work_3/Predict_All_CNVR/model_test_numSNP/functions_gatk/res/res_summary"
res_pred_simple <- readRDS(file = file.path(path_output, "dat_cnvr_all_pred.rds"))
n_sample <- length(unique(res_pred_simple$Sample_ID))
n_cnvr   <- length(unique(res_pred_simple$CNVR_ID))

cutoffs_gq <- c(0, 5, 10, 15, 20, 30, 50) ## change

for (i in 1:length(cutoffs_gq)) {
  
  gq1 <- cutoffs_gq[i]  # cutoff of GQ score
  
  res_pred_simple_gq <- res_pred_simple
  idxs_nocall <- which(res_pred_simple_gq$value_GQ <= gq1)
  if (length(idxs_nocall) >= 1) {
    res_pred_simple_gq$CN_gatk_pred[idxs_nocall] <- 4 ## nocall set as 4
  }
  
  res_pred_simple_gq <- res_pred_simple_gq[order(res_pred_simple_gq$Sample_ID, res_pred_simple_gq$CNVR_ID), ]
  
  m1 <- matrix(res_pred_simple_gq$CN_gatk_pred, nrow = n_cnvr, ncol = n_sample)
  
  freqs <- unlist(lapply(1:nrow(m1), FUN = function(k) {
    v1 <- as.vector(m1[k, ])
    n1 <- length(which(v1 %in% c(0, 1, 3)))
    n1
  }))
  
  
}