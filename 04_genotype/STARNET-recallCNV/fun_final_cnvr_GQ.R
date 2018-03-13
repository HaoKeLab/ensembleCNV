
# final dt_cnvr
path_input <- "/sc/orga/projects/haok01a/chengh04/StarNet/work_3/Predict_All_CNVR/model_test_numSNP/functions_gatk/res/res_summary"
res_pred_simple_raw <- readRDS(file = file.path(path_input, "dat_cnvr_all_pred.rds"))

cutoff_GQ <- 20

idxs_nocall <- which(res_pred_simple_raw$value_GQ <= 20)
res_pred_simple_raw$CN_gatk_pred[idxs_nocall] <- 4  # 4
 
res_pred_simple_raw <- res_pred_simple_raw[order(res_pred_simple_raw$Sample_ID, res_pred_simple_raw$CNVR_ID), ]
n_sample = 834
n_cnvr = 4006
cnvrs <- res_pred_simple_raw$CNVR_ID[1:n_cnvr]
samples <- unique(res_pred_simple_raw$Sample_ID)

mat1 <- matrix(res_pred_simple_raw$CN_gatk_pred, nrow = n_cnvr, ncol = n_sample, dimnames = list(cnvrs, samples))

# cnvrs <- res_pred_simple_raw$CNVR_ID[1:n_cnvr] # all cnvrs

freqs <- unlist(lapply(1:n_cnvr, FUN = function(k) {
  v1 <- as.vector(mat1[k, ])
  sum(v1 %in% c(0, 1, 3))
}))

call_rates <- unlist(lapply(1:n_cnvr, FUN = function(k) {
  v1 <- as.vector(mat1[k, ])
  sum(v1 %in% c(0, 1, 2, 3))/834
}))

# CNVR summary
dat <- data.frame(CNVR_ID = cnvrs, 
                  freq = freqs,
                  call_rate = call_rates, stringsAsFactors = FALSE)

path_output <- "/sc/orga/projects/haok01a/chengh04/StarNet/work_3/Predict_All_CNVR/model_test_numSNP/functions_gatk/res/res_summary/GQ_20_callrate_0"
## filter freq = 0
dat <- subset(dat, freq > 0)
nrow(dat)  ## 2729
sum(dat$freq >= 40)
# save dat_cnvr_info.rds
saveRDS(dat, file = file.path(path_output, "dat_cnvr_info.rds"))

cnvrs1 <- dat$CNVR_ID
idxs_sub <- which(cnvrs %in% cnvrs1)
mat2 <- mat1[idxs_sub, ]
saveRDS(mat2, file = file.path(path_output, "matrix_final.rds"))


dat <- subset(res_pred_simple_raw, CNVR_ID %in% cnvrs1)
nrow(dat)/834
saveRDS(dat, file = file.path(path_output, "eCNV_final.rds"))

file_cnvr <- "cnvrs_boundery_with_batch_info.rds"
path_cnvr <- "/sc/orga/projects/haok01a/chengh04/StarNet/work_3/Predict_All_CNVR/model_test_numSNP/data"
dt_cnvr <- readRDS(file = file.path(path_cnvr, file_cnvr))

dt_cnvr1 <- subset(dt_cnvr, CNVR_ID %in% cnvrs1)
saveRDS(dt_cnvr1, file = file.path(path_output, "eCNV_boundary_final.rds"))


# # set call_rate >= 0.9
# dat1 <- subset(dat, call_rate >= 0.9)  # 3479
# nrow(dat1)
# sum(dat1$freq >= 40) # 138
# 
# cnvrs1 <- dat1$CNVR_ID
# idxs_sub <- which(cnvrs %in% cnvrs1)
# 
# mat2 <- mat1[idxs_sub, ]
# saveRDS(mat2, file = file.path(path_input, "matrix_final_for_compare_dup.rds"))
# 
# file_cnvr <- "cnvrs_boundery_with_batch_info.rds"
# path_cnvr <- "/sc/orga/projects/haok01a/chengh04/StarNet/work_3/Predict_All_CNVR/model_test_numSNP/data"
# dt_cnvr <- readRDS(file = file.path(path_cnvr, file_cnvr))
# 
# dt_cnvr1 <- subset(dt_cnvr, CNVR_ID %in% cnvrs1)
# saveRDS(dt_cnvr1, file = file.path(path_input, "eCNV_boundary_final.rds"))
# 
# saveRDS(dat1, file = file.path(path_input, "cnvrs_info_final.rds"))


