
# add GQ_score as parameter

# transform model to data.frame
trans_model <- function(model) {
  
  mus = model$mu
  sigmas = model$sigma
  lambdas = model$lambda
  
  res <- data.frame(mu0 = mus[1], mu1 = mus[2], mu2 = mus[3], mu3 = mus[4],
                    sigma0 = sigmas[1], sigma1 = sigmas[2], sigma2 = sigmas[3], sigma3 = sigmas[4],
                    lambda0 = lambdas[1], lambda1 = lambdas[2], lambda2 = lambdas[3], lambda3 = lambdas[4])
  
}

# main pipeline function
pipeline_main <- function(dt_cnvrs, paras_LRR, dup_pairs, samples_LRR,
                          plot_steps = TRUE, path_png, GQ_score = 0) {
  
  path_png_diag <- file.path(path_png, "diag")
  path_png_steps <- file.path(path_png, "steps")
  path_png_summary <- file.path(path_png, "summary")
  
  cnvr_id <- unique(dt_cnvrs$CNVR_ID)
  numsnp <- unique(dt_cnvrs$numSNP) # numsnp 
  dt_cnvr <- process_cnvr_LRR(dt_cnvrs = dt_cnvrs, samples_LRR = samples_LRR)
  
  # set CN = 0 cutoff = -0.8
  dt_cnvr0 <- subset(dt_cnvr, LRR_median <= -0.8)
  n0 <- nrow(dt_cnvr0)  # set cutoff of n0 is 5
  
  dt_cnvr_train <- subset(dt_cnvr, LRR_median > - 0.8 & CN != 0) # all confrimed CN = 1/2/3
  
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
  png(filename = file.path(path_png_diag, file_diagnosis), width = 12, height = 12, units = "in", res = 512)
  plot_gmm_diagnosis(dt_cnvr = dt_cnvr_train, paras_model = paras_model)
  dev.off()
  
  # set CN = 0 
  mu0 <- -3
  # sigma0 <- sigma1*10  ## MUST BE CHANGED
  sigma0 <- 0.8*0.8
  if (n0 != 0) {
    if (n0  >= 5) {
      mu0 <- median(dt_cnvr0$LRR_median)
      sigma0 <- sd(dt_cnvr0$LRR_median)
    } 
  }
  
  # split CN = 0 
  # number = 0
  if (n0 <= 5) {
    # model1 <- normalmixEM(x = dt_cnvr$LRR_median, k = 3, 
    #                       mean.constr = c(mu1, mu2, mu3),
    #                       sd.constr = c(sigma1, sigma2, sigma3))
    # add CN = 0 parameters
    model1 <- list()
    model1$mu <- c(mu0, paras_all$mus)
    model1$sigma <- c(sigma0, paras_all$sigmas)
    
    pers <- c(n0/834, ((834-n0)/834)*paras_all$lambdas)
    model1$lambda <- pers
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
  png(filename = file.path(path_png_steps, file_steps), width = 12, height = 12, units = "in", res = 512)
  plot_steps(dt_cnvr_train = dt_cnvr_train, dup_pairs = dup_pairs, dt_cnvr_raw = dt_cnvr, 
             paras = paras_all, dt_LRRBAF = res_gatk_pred1)  ## here
  dev.off()
  
  # if number of CN = 1 <= 2* CN= 0
  # re_model
  n0_new <- sum(res_gatk_pred1$CN_gatk_pred == 0)
  n1_new <- sum(res_gatk_pred1$CN_gatk_pred == 1)
  # n0_new
  # n1_new
  
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
  png(filename = file.path(path_png_steps, file_steps), width = 12, height = 12, units = "in", res = 512)
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
  
  # add GQ_score cutoff here
  idxs_nocall <- which(res_gatk_pred_final$value_GQ <= GQ_score)
  call_rate <- 1 - length(idxs_nocall)/nrow(res_gatk_pred_final)
  if (length(idxs_nocall) >= 1) {
    res_gatk_pred_final$CN_gatk_pred[idxs_nocall] <- 4
  }
  
  # scatter plot LRR_median
  dt_pred <- res_gatk_pred_final[, c("Sample_ID", "CN_gatk_pred")]
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
  
  # plot_raw
  # add gray color point here
  plot_gatk <- ggplot() + 
    geom_point(data = subset(dt_cnvr_scatter, CN_gatk_pred == 0), aes(idx, LRR_median), col = "black") +
    geom_point(data = subset(dt_cnvr_scatter, CN_gatk_pred == 1), aes(idx, LRR_median), col = "red") + 
    geom_point(data = subset(dt_cnvr_scatter, CN_gatk_pred == 2), aes(idx, LRR_median), col = "green") +
    geom_point(data = subset(dt_cnvr_scatter, CN_gatk_pred == 3), aes(idx, LRR_median), col = "blue") +
    geom_point(data = subset(dt_cnvr_scatter, CN_gatk_pred == 4), aes(idx, LRR_median), col = "gray") +
    theme_bw(base_size = 10) +
    ggtitle(label = "CNV call from gatk similary method",
            subtitle = paste("GQ score:", GQ_score, "call rate:", round(call_rate, 3)))
  
  # plot_gatk
  
  filefinal <- paste0("summary_", cnvr_id, ".png")
  png(filename = file.path(path_png_summary, filefinal),
      width = 12, height = 12, units = "in", res = 512)
  grid.arrange(plot_MAF, plot_final, plot_raw, plot_gatk, nrow = 2)
  dev.off()

  res_pars <- trans_model(model = model_final)
  res_pars$CNVR_ID = cnvr_id
  res_pars$numSNP  = numsnp
  
  # res_gatk_pred_final  ## final result
  res <- list(res_gatk_pred_final = res_gatk_pred_final,
              res_pars = res_pars)  ## paras for each CNVR_ID
}
 


