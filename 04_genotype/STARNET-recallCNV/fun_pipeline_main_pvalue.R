
rm(list = ls())

# set cutoff of GQ-score 
# 20/50/
# result save png_pvalue
setwd("C:/mssm_projects/STARNET/work_3/Predict_All_CNVR/model_test_numSNP/functions_gatk/")

# load packages
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

source('./fun_BAF.R')
source('./fun_gatk.R')
source('./fun_LRR.R')
source('./fun_models_zz_8_cc_raw.R')
source('./fun_par.R')
source('./fun_plot_steps.R')
source('./fun_plot_diagnosis.R')
source('./fun_heatmap.R')

# main pipeline function
pipeline_main <- function(dt_cnvrs, paras_LRR, dup_pairs, samples_LRR,
                          plot_steps = TRUE, cutoff_GQ_score) {
  
  cnvr_id <- unique(dt_cnvrs$CNVR_ID)
  numsnp <- unique(dt_cnvrs$numSNP) # numsnp 
  dt_cnvr <- process_cnvr_LRR(dt_cnvrs = dt_cnvrs, samples_LRR = samples_LRR)
  
  # set CN = 0 cutoff = -0.8
  dt_cnvr0 <- subset(dt_cnvr, LRR_median <= -0.8)
  n0 <- nrow(dt_cnvr0)  # set cutoff of n0 is 5
  
  dt_cnvr_train <- subset(dt_cnvr, LRR_median > - 0.8) # all confrimed CN = 1/2/3
  
  res_paras <- train_model_zz(dt_cnvr = dt_cnvr_train, paras_LRR = paras_LRR)

  paras_all <- res_paras$paras_all  # final paras from gmm model
  paras_model <- res_paras$paras_model  # step paras for plot diagnosis
  
  ## all predict result
  mu1 <- paras_all$mus[1]
  sigma1 <- paras_all$sigmas[1]
  mu2 <- paras_all$mus[2]
  sigma2 <- paras_all$sigmas[2]
  mu3 <- paras_all$mus[3]
  sigma3 <- paras_all$sigmas[3]
  
  # save diagnosis png
  ## plot diagnosis 
  file_diagnosis <- paste0("diag_", cnvr_id, ".png")
  png(filename = paste0('./png_pvalue/', file_diagnosis), width = 12, height = 12, units = "in", res = 512)
  plot_gmm_diagnosis(dt_cnvr = dt_cnvr_train, paras_model = paras_model)
  dev.off()
  
  # set CN = 0 
  mu0 <- -3
  sigma0 <- sigma1*10  ## MUST BE CHANGED
  if (n0 != 0) {
    mu0 <- median(dt_cnvr0$LRR_median)
    if (n0  >= 5) {
      sigma0 <- sd(dt_cnvr0$LRR_median)
    }
  }
  
  # split CN = 0 
  # number = 0
  if (n0 == 0) {
    model1 <- normalmixEM(x = dt_cnvr$LRR_median, k = 3, 
                          mean.constr = c(mu1, mu2, mu3),
                          sd.constr = c(sigma1, sigma2, sigma3))
    # add CN = 0 parameters
    model1$mu <- c(mu0, model1$mu)
    model1$sigma <- c(sigma0, model1$sigma)
    model1$lambda <- c(0, model1$lambda)
  } else {
    # number of CN = 1
    model1 <- normalmixEM(x = dt_cnvr$LRR_median, k = 4,
                          mean.constr = c(mu0, mu1, mu2, mu3),
                          sd.constr = c(sigma0, sigma1, sigma2, sigma3))
  }
  
  # calculate BAF
  dt_BAF1 <- calculate_BAF_gatk_whole(dt_cnvrs = dt_cnvrs)
  # calculate LRR
  dt_LRR1 <- output_LRR_gatk(dt_cnvr = dt_cnvr, model = model1)
  dt_LRRBAF1 <- merge(dt_LRR1, dt_BAF1)
  res_gatk_pred1 <- output_gatk_result(dt_LRRBAF = dt_LRRBAF1)
  mean1_GQ <- mean(res_gatk_pred1$value_GQ)
  
  ## save steps_1_ png
  file_steps <- paste0("steps_1_", cnvr_id, ".png")
  png(filename = paste0('./png_pvalue/', file_steps), width = 12, height = 12, units = "in", res = 512)
  plot_steps(dt_cnvr_train = dt_cnvr_train, dup_pairs = dup_pairs, dt_cnvr_raw = dt_cnvr, 
             paras = paras_all, dt_LRRBAF = res_gatk_pred1)  ## here
  dev.off()
  
  # if number of CN = 1 <= 2* CN= 0
  # re_model
  n0_new <- sum(res_gatk_pred1$CN_gatk_pred == 0)
  n1_new <- sum(res_gatk_pred1$CN_gatk_pred == 1)
  n0_new
  n1_new
  
  # hardy weinberg test
  res_gatk_pred_final <- NULL
  model_final <- NULL
  if (n1_new >= n0_new) {
    res_gatk_pred_final <- res_gatk_pred1
    model_final <- model1  ## final model
  } else {
    mu1 <- paras_all$mus[2]
    sigma1 <- paras_all$sigmas[2]
    mu2 <- paras_all$mus[3]
    sigma2 <- paras_all$sigmas[3]
    mu3 <- paras_all$mus[3] + paras_LRR$LRR_mean$CN_3
    sigma3 <- paras_LRR$LRR_sd$CN_3/sqrt(numsnp)
    
    if (n0 == 0) {
      model2 <- normalmixEM(x = dt_cnvr$LRR_median, k = 3, 
                            mean.constr = c(mu1, mu2, mu3),
                            sd.constr = c(sigma1, sigma2, sigma3))
      # add CN = 0 parameters
      model2$mu <- c(mu0, model2$mu)
      model2$sigma <- c(sigma0, model2$sigma)
      model2$lambda <- c(0, model2$lambda)
    } else {
      model2 <- normalmixEM(x = dt_cnvr$LRR_median, k = 4,
                            mean.constr = c(mu0, mu1, mu2, mu3),
                            sd.constr = c(sigma0, sigma1, sigma2, sigma3))
    }
    
    
    # calculate BAF
    dt_BAF2 <- calculate_BAF_gatk_whole(dt_cnvrs = dt_cnvrs)
    # calculate LRR
    dt_LRR2 <- output_LRR_gatk(dt_cnvr = dt_cnvr, model = model2)
    dt_LRRBAF2 <- merge(dt_LRR2, dt_BAF2)
    res_gatk_pred2 <- output_gatk_result(dt_LRRBAF = dt_LRRBAF2)
    mean2_GQ <- mean(res_gatk_pred2$value_GQ)
    
    cat(mean1_GQ, mean2_GQ, "\n")
    
    if (mean2_GQ > mean1_GQ) {
      res_gatk_pred_final <- res_gatk_pred2
      model_final <- model2  ## model final
    } else {
      res_gatk_pred_final <- res_gatk_pred1
      model_final <- model1  ## model final
    }
  }
  
  ## save steps_2 png
  file_steps <- paste0("steps_2_", cnvr_id, ".png")
  png(filename = paste0('./png_pvalue/', file_steps), width = 12, height = 12, units = "in", res = 512)
  plot_steps(dt_cnvr_train = dt_cnvr_train, dup_pairs = dup_pairs, dt_cnvr_raw = dt_cnvr, 
             paras = paras_all, dt_LRRBAF = res_gatk_pred_final)  ## here
  dev.off()
  
  # plot for final
  # for new input must ordered as follow
  # dt_pfb <- dt_cnvrs[order(dt_cnvrs$Sample_ID, dt_cnvrs$Position), ]
  dt_pfb <- dt_cnvrs[1:numsnp, ]
  dt_pfb$MAF <- pmin(dt_pfb$PFB, 1 - dt_pfb$PFB) ##
  dt_pfb <- dt_pfb[, c("Name", "MAF")]
  
  plot_MAF <- ggplot(data = dt_pfb, aes(Name, MAF)) + 
    geom_col() + 
    ggtitle(label = paste("snps MAF in Position order", "numSNP:", numsnp)) + 
    labs(x = "SNP Name") + 
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # final 
  plot_final <- plot_model_final(paras = model_final, dt_cnvr = dt_cnvr_train,
                                 title = paste("final model for", cnvr_id, "numSNP:", numsnp))
  
  # scatter plot LRR_median
  dt_pred <- res_gatk_pred_final[, c("Sample_ID", "CN_gatk_pred", "value_GQ")]
  
  # set value_GQ
  idx_nocall <- which(dt_pred$value_GQ <= cutoff_GQ_score)
  if (length(idx_nocall) >= 1) {
    dt_pred$CN_gatk_pred[idx_nocall] <- 4
  }
  
  dt_cnvr_scatter <- merge(dt_cnvr, dt_pred)
  dt_cnvr_scatter <- dt_cnvr_scatter[order(dt_cnvr_scatter$CN), ]
  dt_cnvr_scatter$idx <- 1:nrow(dt_cnvr_scatter)
  myColors <- brewer.pal(4, "Set1")
  plot_raw <- ggplot() + 
    geom_point(data = subset(dt_cnvr_scatter, CN == 0), aes(idx, LRR_median), col = "black") +
    geom_point(data = subset(dt_cnvr_scatter, CN == 1), aes(idx, LRR_median), col = "red") + 
    geom_point(data = subset(dt_cnvr_scatter, CN == 2), aes(idx, LRR_median), col = "green") +
    geom_point(data = subset(dt_cnvr_scatter, CN == 3), aes(idx, LRR_median), col = "blue") +
    theme_bw(base_size = 10) + 
    ggtitle(label = "CNV call from IPQ") 
  
  plot_raw
  
  plot_gatk <- ggplot() + 
    geom_point(data = subset(dt_cnvr_scatter, CN_gatk_pred == 0), aes(idx, LRR_median), col = "black") +
    geom_point(data = subset(dt_cnvr_scatter, CN_gatk_pred == 1), aes(idx, LRR_median), col = "red") + 
    geom_point(data = subset(dt_cnvr_scatter, CN_gatk_pred == 2), aes(idx, LRR_median), col = "green") +
    geom_point(data = subset(dt_cnvr_scatter, CN_gatk_pred == 3), aes(idx, LRR_median), col = "blue") +
    geom_point(data = subset(dt_cnvr_scatter, CN_gatk_pred == 4), aes(idx, LRR_median), col = "gray") +
    theme_bw(base_size = 10) +
    ggtitle(label = "CNV call from gatk similary method") 
  
  plot_gatk
  
  filefinal <- paste0("summary_", cnvr_id, ".png")
  png(filename = file.path('./png_pvalue/', filefinal),
      width = 12, height = 12, units = "in", res = 512)
  grid.arrange(plot_MAF, plot_final, plot_raw, plot_gatk, nrow = 2)
  dev.off()

  res_gatk_pred_final  ## final result
}
 

# path_input <- "C:/mssm_projects/STARNET/work_3/Predict_All_CNVR/model_test_numSNP/res_CNVR_with_CN1/dat_CNVR_check/"
# files <- list.files(path = path_input)
# 
# for (f1 in files) {
#   
#   cat(f1, "\n")
#   filename <- file.path(path_input, f1)
#   dt_cnvrs <- readRDS(file = filename)
#   
#   res <- pipeline_main(dt_cnvrs = dt_cnvrs, paras_LRR = paras_LRR,
#                        dup_pairs = dup_pairs, samples_LRR = samples_LRR,
#                        plot_steps = TRUE)
# }
# 





filename <- "C:/mssm_projects/STARNET/work_3/Predict_All_CNVR/model_test_numSNP/res_CNVR_with_CN1/dat_CNVR_check/CNVR_79_r1_chr1_p.rds "
dt_cnvrs <- readRDS(file = filename)

table(dt_cnvrs$CN)/unique(dt_cnvrs$numSNP)

res <- pipeline_main(dt_cnvrs = dt_cnvrs, paras_LRR = paras_LRR,
                     dup_pairs = dup_pairs, samples_LRR = samples_LRR,
                     plot_steps = TRUE, cutoff_GQ_score = 20)

hist(res$value_GQ, breaks = 50)
mean(res$value_GQ) # 191.4534
table(res$CN, res$CN_gatk_pred)
