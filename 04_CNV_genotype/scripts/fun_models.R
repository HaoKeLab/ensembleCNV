
mode_zz <- function(dt_cnvr) {
  dt <- dt_cnvr
  cutoffs <- quantile(dt$LRR_median, c(0.1, 0.9))
  dt1 <- dt[which(dt$LRR_median >= cutoffs[1] & dt$LRR_median <= cutoffs[2]), ]
  md1 <- mlv(x = dt1$LRR_median, method = "parzen", kernel = "gaussian")
  #m1  <- md1$M
  
  sd1 <- sd(dt1$LRR_median)
  
  return(list(mu = md1, sigma = sd1))
}

# CN = 1 and 3 all have more 10 samples
model1_zz <- function(dt_cnvr, paras_LRR) {
  
  numsnp <- unique(dt_cnvr$numSNP)
  # init paras_model
  paras_model <- list()
  
  paras_model$stage1 <- list()
  paras_model$stage2 <- list()
  
  dt_cnvr <- dt_cnvr[order(dt_cnvr$CN), ]
  n <- nrow(dt_cnvr)  # number of samples
  numsnp <- unique(dt_cnvr$numSNP) # numsnp
  cnvr_id <- unique(dt_cnvr$CNVR_ID)
  chr <- unique(dt_cnvr$Chr)
  
  tbl1 <- table(dt_cnvr$CN)
  # first round gmm model
  lambdas1 <- prop.table(tbl1)
  
  model1 <- list()
  
  n1s <- as.vector(tbl1)  ##
  # add here percent of CN = 2 < 50%
  if (lambdas1[2] <= 0.5) {
    
    # exclude P/Q/PQ
    idxs <- which(dt_cnvr$alg %in% c("I", "IP", "IPQ", "IQ"))
    
    dt_new <- dt_cnvr[idxs, ]
    dt2_new <- NULL
    if (length(idxs) == 0) {
      dt2_new <- dt_cnvr
    } else {
      dt2_new <- dt_cnvr[-idxs, ]
    }
    
    
    m2 <- mode_zz(dt_cnvr = dt2_new)
    mu2 <- m2$mu
    sigma2 <- m2$sigma
    
    dt1_new <- subset(dt_new, CN == 1)
    dt3_new <- subset(dt_cnvr, CN == 3)
    
    n1s <- c(nrow(dt1_new), nrow(dt2_new), nrow(dt3_new)) ## n1s
    
    if (nrow(dt1_new) >= 60) {
      mu1 <- mean(dt1_new$LRR_median)
      sigma1 <- sd(dt1_new$LRR_median)
    } else {
      mu1 <- mu2 + paras_LRR$LRR_mean$CN_1
      sigma1 <- paras_LRR$LRR_sd$CN_1/sqrt(numsnp)
    }
    
    if (nrow(dt3_new) >= 60) {
      mu3 <- mean(dt3_new$LRR_median)
      sigma3 <- sd(dt3_new$LRR_median)
    } else {
      mu3 <- mu2 + paras_LRR$LRR_mean$CN_3
      sigma3 <- paras_LRR$LRR_sd$CN_3/sqrt(numsnp)
    } 
    
    if ((mu2 - mu1) < abs(0.2*paras_LRR$LRR_mean$CN_1)) {
      
      mu1 <- mu2 + paras_LRR$LRR_mean$CN_1
      
      num_permute <- 60
      
      sigma1 <- paras_LRR$LRR_sd$CN_1/sqrt(numsnp)
      LRR_median_cn1_permute <- rnorm(num_permute, mean = mu1, sd = sigma1)
      
      dt_cn1_permute <- data.frame(Sample_ID = paste0("Sample_CN1_", 1:num_permute),
                                   CNVR_ID = rep(cnvr_id, num_permute),
                                   LRR_median = LRR_median_cn1_permute,
                                   Chr = chr, alg = rep("permute", num_permute),
                                   CN = 1, numSNP = numsnp, stringsAsFactors = FALSE)
      
      dt_cnvr <- rbind(dt_cnvr, dt_cn1_permute)
      
      n1s[1] <- n1s[1] + 60
    }
    
    lambdas1 <- n1s/sum(n1s)
    
    mus1 <- c(mu1, mu2, mu3)
    sigmas1 <- c(sigma1, sigma2, sigma3)
    # init parameters
    paras_model$stage1$init <- list(mu = c(mu1, mu2, mu3),
                                    sigma = c(sigma1, sigma2, sigma3),
                                    lambda = lambdas1)
    
    if (nrow(dt1_new) < 10 & nrow(dt3_new) < 10) {
      model1 <- normalmixEM(x = dt_cnvr$LRR_median, mu = mus1, 
                            mean.constr = c(mu1, NA, mu3), 
                            sigma = sigmas1, k = 3)
    } else if (nrow(dt1_new) < 10 & nrow(dt3_new) >= 10) {
      model1 <- normalmixEM(x = dt_cnvr$LRR_median, mu = mus1, 
                            mean.constr = c(mu1, NA, NA),
                            sigma = sigmas1, k = 3)
    } else if (nrow(dt1_new >= 10) & nrow(dt3_new) < 10) {
      model1 <- normalmixEM(x = dt_cnvr$LRR_median, mu = mus1, 
                            mean.constr = c(NA, NA, mu3),
                            sigma = sigmas1, k = 3)
    } else {
      model1 <- normalmixEM(x = dt_cnvr$LRR_median, mu = mus1, 
                            sigma = sigmas1, k = 3)
    }
    
  } else  {
    
    dt2 <- subset(dt_cnvr, CN == 2)
    paras2 <- mode_zz(dt_cnvr = dt2)
    mu2 <- paras2$mu
    sigma2 <- paras2$sigma
    # CN = 1
    dt1 <- subset(dt_cnvr, CN == 1)
    mu1 <- mean(dt1$LRR_median)
    sigma1 <- sd(dt1$LRR_median)
    
    # CN = 3
    dt3 <- subset(dt_cnvr, CN == 3)
    mu3 <- mean(dt3$LRR_median)
    sigma3 <- sd(dt3$LRR_median)
    
    if ((mu2 - mu1) < abs(0.2*paras_LRR$LRR_mean$CN_1)) {
      mu1 <- mu2 + paras_LRR$LRR_mean$CN_1
      
      num_permute <- 60
      
      sigma1 <- paras_LRR$LRR_sd$CN_1/sqrt(numsnp)
      LRR_median_cn1_permute <- rnorm(num_permute, mean = mu1, sd = sigma1)
      
      dt_cn1_permute <- data.frame(Sample_ID = paste0("Sample_CN1_", 1:num_permute),
                                   CNVR_ID = rep(cnvr_id, num_permute),
                                   LRR_median = LRR_median_cn1_permute,
                                   Chr = chr, alg = rep("permute", num_permute),
                                   CN = 1, numSNP = numsnp, stringsAsFactors = FALSE)
      
      dt_cnvr <- rbind(dt_cnvr, dt_cn1_permute)
      
      n1s[1] <- n1s[1] + num_permute
    }
    
    mus1 <- c(mu1, mu2, mu3)
    sigmas1 <- c(sigma1, sigma2, sigma3)
    
    lambdas1 <- n1s/sum(n1s)
    # init parameters
    paras_model$stage1$init <- list(mu = c(mu1, mu2, mu3),
                                    sigma = c(sigma1, sigma2, sigma3),
                                    lambda = lambdas1)
    
    model1 <- normalmixEM(x = dt_cnvr$LRR_median, mu = mus1, 
                          sigma = sigmas1, lambda = lambdas1, k = 3)
    
  }
  
  
  # model parameters
  paras_model$stage1$model <- list(mu = model1$mu,
                                   sigma = model1$sigma,
                                   lambda = model1$lambda)
  
  #--------------------------------------------
  cns_pred <- apply(model1$posterior, MARGIN = 1, which.max)
  probs <- model1$posterior
  # CN = 1
  probs1 <- as.vector(probs[, 1]) 
  idxs1  <- which(cns_pred == 1 & probs1 >= 0.9)
  n1 <- length(idxs1)
  n1
  # CN = 3
  probs3 <- as.vector(probs[, 3])
  idxs3  <- which(cns_pred == 3 & probs3 >= 0.9)
  n3 <- length(idxs3)
  n3
  # CN = 2
  idxs2 <- setdiff(1:nrow(dt_cnvr), union(idxs1, idxs3))
  dt22 <- dt_cnvr[idxs2, ]
  m2 <- mode_zz(dt_cnvr = dt22)
  mu2 <- m2$mu
  sigma2 <- m2$sigma
  n2 <- nrow(dt22)
  
  # add 
  n2s <- c(n1, n2, n3)
  flag1 <- ifelse(n2s[1]/n1s[1] >= 0.8,TRUE,FALSE)
  flag3 <- ifelse(n2s[3]/n1s[3] >= 0.8,TRUE,FALSE)  ## prop
  
  # parameters
  mu1 <- NA
  sigma1 <- NA
  percent.dist12 = abs(model1$mu[2] - model1$mu[1])/abs(paras_LRR$LRR_mean$CN_1)
  if ( percent.dist12 <= 0.2) {
    mu1 <- mu2 + paras_LRR$LRR_mean$CN_1
    sigma1 <- paras_LRR$LRR_sd$CN_1/sqrt(numsnp)
    n1 = 0
  } else {
    if (n1 >= 10) {
       mu1 <- mean(dt_cnvr$LRR_median[idxs1])
       sigma1 <- sd(dt_cnvr$LRR_median[idxs1])
    } else {
       mu1 <- mu2 + paras_LRR$LRR_mean$CN_1
       sigma1 <- paras_LRR$LRR_sd$CN_1/sqrt(numsnp)
    }
  }
  
  
  
  mu3 <- NA
  sigma3 <- NA
  if (n3 >= 10) {
    mu3 <- mean(dt_cnvr$LRR_median[idxs3])
    sigma3 <- sd(dt_cnvr$LRR_median[idxs3])
  } else {
    mu3 <- mu2 + paras_LRR$LRR_mean$CN_3
    sigma3 <- paras_LRR$LRR_sd$CN_3/sqrt(numsnp)
  }
  
  if (mu3 < mu2) {
    mu3 <- mu2 + paras_LRR$LRR_mean$CN_3
    sigma3 <- paras_LRR$LRR_sd$CN_3/sqrt(numsnp)
  }
  
  lambdas2 <- c(n1, n2, n3)/n
  paras_model$stage2$init <- list(mu = c(mu1, mu2, mu3),
                                  sigma = c(sigma1, sigma2, sigma3),
                                  lambda = lambdas2)
  model2 <- list()
  if (n1 >= 10 & n3 >= 10) {
    model2 <- normalmixEM(x = dt_cnvr$LRR_median, mu = c(mu1, mu2, mu3),
                          sigma = c(sigma1, sigma2, sigma3), lambda = lambdas2,
                          k = 3)
  } else if (n1 >= 10 & n3 < 10) {
    model2 <- normalmixEM(x = dt_cnvr$LRR_median, mu = c(mu1, mu2, mu3),
                          sigma = c(sigma1, sigma2, sigma3), 
                          mean.constr = c(NA, NA, mu3), k = 3)
  } else if (n1 < 10 & n3 >= 10) {
    model2 <- normalmixEM(x = dt_cnvr$LRR_median, mu = c(mu1, mu2, mu3),
                          sigma = c(sigma1, sigma2, sigma3), 
                          mean.constr = c(mu1, NA, mu3), k = 3)
  } else if (n1 < 10 & n3 < 10) {
    mu2 <- median(dt22$LRR_median)
    sigma2 <- sd(dt22$LRR_median)
    
    mu1 <- mu2 + paras_LRR$LRR_mean$CN_1
    if(n1 == 0) {
      sigma1 <- Inf
    } else {
      sigma1 <- paras_LRR$LRR_sd$CN_1/sqrt(numsnp)
    }
    
    mu3 <- mu2 + paras_LRR$LRR_mean$CN_3
    if (n3 == 0) {
      sigma3 <- Inf
    } else {
      sigma3 <- paras_LRR$LRR_sd$CN_3/sqrt(numsnp)
    }
    
    model2$mu <- c(mu1, mu2, mu3)
    model2$sigma <- c(sigma1, sigma2, sigma3)
    model2$lambda <- c(1, 1, 1)
  }
  
  # check parameters 
  cutoff1 <- 0.5
  cutoff3 <- 0.8
  mu1f <- model2$mu[1]
  mu2f <- model2$mu[2]
  mu3f <- model2$mu[3]
  
  f1 <- (abs(mu1f - mu2f) > cutoff1*abs(paras_LRR$LRR_mean$CN_1))
  f3 <- (abs(mu2f - mu3f) > cutoff3*abs(paras_LRR$LRR_mean$CN_3))
  
  f1 <- f1|flag1
  f3 <- f3|flag3
  
  if (f1 == FALSE & f3 == FALSE) {
    
    mu2 <- median(dt_cnvr$LRR_median)
    sd2 <- sd(dt_cnvr$LRR_median)
    
    mu3 <- mu2 + paras_LRR$LRR_mean$CN_3
    mu1 <- mu2 + paras_LRR$LRR_mean$CN_1
    
    model3 <- normalmixEM(x = dt_cnvr$LRR_median, mu = c(mu1, mu2, mu3), 
                          mean.constr = c(mu1, NA, mu3), k = 3)
    
    paras_model$stage2$model <- list(mu = model3$mu,
                                     sigma = model3$sigma,
                                     lambda = model3$lambda)
    
    paras_all <- list(mus = model3$mu, 
                      sigmas = model3$sigma,
                      lambdas = model3$lambda)
    
    
  } else if (f1 == FALSE & f3 == TRUE) {
    
    mu1 <- mu2 + paras_LRR$LRR_mean$CN_1
    
    model3 <- normalmixEM(x = dt_cnvr$LRR_median, mu = c(mu1, mu2, mu3),
                          mean.constr = c(mu1, NA, mu3), k = 3)
    
    paras_model$stage2$model <- list(mu = model3$mu,
                                     sigma = model3$sigma,
                                     lambda = model3$lambda)
    
    paras_all <- list(mus = model3$mu, 
                      sigmas = model3$sigma,
                      lambdas = model3$lambda)
    
  } else if (f1 == TRUE & f3 == FALSE) {
    
    mu3 <- mu2 + paras_LRR$LRR_mean$CN_3
    
    model3 <- normalmixEM(x = dt_cnvr$LRR_median, mu = c(mu1, mu2, mu3),
                          mean.constr = c(mu1, NA, mu3), k = 3)
    
    paras_model$stage2$model <- list(mu = model3$mu,
                                     sigma = model3$sigma,
                                     lambda = model3$lambda)
    
    paras_all <- list(mus = model3$mu, 
                      sigmas = model3$sigma,
                      lambdas = model3$lambda)
    
  } else {
    
    mu22 <- NA
    if (n1 <= 5 & n3 <= 5) {
      mu22 <- mu2
    }
    model3 <- normalmixEM(x = dt_cnvr$LRR_median, mu = c(mu1, mu2, mu3),
                          mean.constr = c(mu1, mu22, mu3), k = 3)
    
    paras_model$stage2$model <- list(mu = model3$mu,
                                     sigma = model3$sigma,
                                     lambda = model3$lambda)
    
    paras_all <- list(mus = model3$mu, 
                      sigmas = model3$sigma,
                      lambdas = model3$lambda)
    
  }
  
  res <- list(paras_all = paras_all,
              paras_model = paras_model)
  return(res)  # return
}

model2_zz <- function(dt_cnvr, paras_LRR) {
  
  # init paras_model
  paras_model <- list()
  
  paras_model$stage1 <- list()
  paras_model$stage2 <- list()
  
  dt_cnvr <- dt_cnvr[order(dt_cnvr$CN), ]
  n <- nrow(dt_cnvr)  # number of samples
  numsnp <- unique(dt_cnvr$numSNP) # numsnp
  cnvr_id <- unique(dt_cnvr$CNVR_ID)
  chr <- unique(dt_cnvr$Chr)
  
  cn_factor <- factor(dt_cnvr$CN, levels = c(1, 2, 3))
  v <- as.vector(table(cn_factor))
  lambdas1 <- v/n ## 

  n1s <- v
  
  dt1 <- subset(dt_cnvr, CN == 1)
  dt2 <- subset(dt_cnvr, CN == 2)
  dt3 <- subset(dt_cnvr, CN == 3)
  
  mu1 <- mean(dt1$LRR_median)
  sigma1 <- sd(dt1$LRR_median)
  
  paras1 <- mode_zz(dt_cnvr = dt2)
  mu2 <- paras1$mu
  sigma2 <- paras1$sigma
  
  mu3 <- mean(dt3$LRR_median)
  sigma3 <- sd(dt3$LRR_median)
  
  # add 8_cc here
  if ((mu2 - mu1) < abs(0.2*paras_LRR$LRR_mean$CN_1)) {
    mu1 <- mu2 + paras_LRR$LRR_mean$CN_1
    
    num_permute <- 20
    # num_permute <- 60
    
    sigma1 <- paras_LRR$LRR_sd$CN_1/sqrt(numsnp)
    LRR_median_cn1_permute <- rnorm(num_permute, mean = mu1, sd = sigma1)
    
    dt_cn1_permute <- data.frame(Sample_ID = paste0("Sample_CN1_", 1:num_permute),
                                 CNVR_ID = rep(cnvr_id, num_permute),
                                 LRR_median = LRR_median_cn1_permute,
                                 Chr = chr, alg = rep("permute", num_permute),
                                 CN = 1, numSNP = numsnp, stringsAsFactors = FALSE)
    
    dt_cnvr <- rbind(dt_cnvr, dt_cn1_permute)

    n1s[1] <- n1s[1] + num_permute
  }
  
  lambdas1 <- n1s/sum(n1s)  ## add
  
  paras_model$stage1$init <- list(mu = c(mu1, mu2, mu3),
                                  sigma = c(sigma1, sigma2, sigma3),
                                  lambda = lambdas1)

  model1 <- normalmixEM(x = dt_cnvr$LRR_median, mu = c(mu1, mu2, mu3),
                        sigma = c(sigma1, sigma2, sigma3), k = 3)
  
  paras_model$stage1$model <- list(mu = model1$mu,
                                   sigma = model1$sigma,
                                   lambda = model1$lambda)
  
  cns_pred <- apply(model1$posterior, MARGIN = 1, which.max)
  probs    <- model1$posterior
  # CN = 1
  probs1   <- as.vector(probs[, 1])
  idxs1    <- which(probs1 >= 0.9 & cns_pred == 1)
  n1 <- length(idxs1)
  n1
  
  # CN = 3
  probs3 <- as.vector(probs[, 3])
  idxs3 <- which(probs3 >= 0.9 & cns_pred == 3)
  n3 <- length(idxs3)
  n3
  
  idxs <- union(idxs1, idxs3)
  idxs2 <- setdiff(1:nrow(dt_cnvr), idxs)
  dt22 <- dt_cnvr[idxs2, ]
  n2 <- nrow(dt22)
  
  paras2 <- mode_zz(dt_cnvr = dt22)
  mu2 <- paras2$mu
  sigma2 <- paras2$sigma
  
  n2s <- c(n1, n2, n3)
  flag1 <- ifelse(n2s[1]/n1s[1] >= 0.8, TRUE, FALSE)
  flag3 <- ifelse(n2s[3]/n1s[3] >= 0.8, TRUE, FALSE)
  
  # parameters
  mu1 <- NA
  sigma1 <- NA
  if (n1 >= 10) {
    mu1 <- mean(dt_cnvr$LRR_median[idxs1])
    sigma1 <- sd(dt_cnvr$LRR_median[idxs1])
  } else {
    mu1 <- mu2 + paras_LRR$LRR_mean$CN_1
    sigma1 <- paras_LRR$LRR_sd$CN_1
  }
  
  mu3 <- NA
  sigma3 <- NA
  if (n3 >= 10) {
    mu3 <- mean(dt_cnvr$LRR_median[idxs3])
    sigma3 <- sd(dt_cnvr$LRR_median[idxs3])
  } else {
    mu3 <- mu2 + paras_LRR$LRR_mean$CN_3
    sigma3 <- paras_LRR$LRR_sd$CN_3
  }
  ## add
  if (mu3 < mu2) {
    mu3 <- mu2 + paras_LRR$LRR_mean$CN_3
    sigma3 <- paras_LRR$LRR_sd$CN_3
  }
  
  n2 <- n - n1 - n3
  lambdas2 <- c(n1, n2, n3)/n
  # init parameters
  paras_model$stage2$init <- list(mu = c(mu1, mu2, mu3),
                                  sigma = c(sigma1, sigma2, sigma3),
                                  lambda = lambdas2)
  
  model2 <- list()
  if (n1 >= 10 & n3 >= 10) {
    model2 <- normalmixEM(x = dt_cnvr$LRR_median, mu = c(mu1, mu2, mu3),
                          sigma = c(sigma1, sigma2, sigma3), 
                          mean.constr = c(mu1, NA, mu3),
                          k = 3)
  } else if (n1 >= 10 & n3 < 10) {
    model2 <- normalmixEM(x = dt_cnvr$LRR_median, mu = c(mu1, mu2, mu3),
                          sigma = c(sigma1, sigma2, sigma3), 
                          mean.constr = c(mu1, NA, mu3), k = 3)
  } else if (n1 < 10 & n3 >= 10) {
    model2 <- normalmixEM(x = dt_cnvr$LRR_median, mu = c(mu1, mu2, mu3),
                          sigma = c(sigma1, sigma2, sigma3), 
                          mean.constr = c(mu1, NA, mu3), k = 3)
  } else if (n1 < 10 & n3 < 10) {
    mu2 <- median(dt22$LRR_median)
    sigma2 <- sd(dt22$LRR_median)
    
    mu1 <- mu2 + paras_LRR$LRR_mean$CN_1
    if(n1 == 0) {
      sigma1 <- Inf
    } else {
      sigma1 <- paras_LRR$LRR_sd$CN_1/sqrt(numsnp)
    }
    
    mu3 <- mu2 + paras_LRR$LRR_mean$CN_3
    if (n3 == 0) {
      sigma3 <- Inf
    } else {
      sigma3 <- paras_LRR$LRR_sd$CN_3/sqrt(numsnp)
    }
    
    model2$mu <- c(mu1, mu2, mu3)
    model2$sigma <- c(sigma1, sigma2, sigma3)
    model2$lambda <- c(1, 1, 1)
  }
  
  
  # check parameters 
  cutoff1 <- 0.5
  cutoff3 <- 0.8
  mu1f <- model2$mu[1]
  mu2f <- model2$mu[2]
  mu3f <- model2$mu[3]
  
  f1 <- (abs(mu1f - mu2f) > cutoff1*abs(paras_LRR$LRR_mean$CN_1))
  f3 <- (abs(mu2f - mu3f) > cutoff3*abs(paras_LRR$LRR_mean$CN_3))
  
  ## 
  f1 <- f1|flag1
  f3 <- f3|flag3
  
  if (f1 == FALSE & f3 == FALSE) {
    
    mu2 <- median(dt_cnvr$LRR_median)
    sd2 <- sd(dt_cnvr$LRR_median)
    
    mu3 <- mu2 + paras_LRR$LRR_mean$CN_3
    mu1 <- mu2 + paras_LRR$LRR_mean$CN_1
    
    model3 <- normalmixEM(x = dt_cnvr$LRR_median, mu = c(mu1, mu2, mu3), 
                          mean.constr = c(mu1, NA, mu3), k = 3)
    
    paras_model$stage2$model <- list(mu = model3$mu,
                                     sigma = model3$sigma,
                                     lambda = model3$lambda)
    
    paras_all <- list(mus = model3$mu, 
                      sigmas = model3$sigma,
                      lambdas = model3$lambda)
    
    
  } else if (f1 == FALSE & f3 == TRUE) {
    
    mu1 <- mu2 + paras_LRR$LRR_mean$CN_1
    
    model3 <- normalmixEM(x = dt_cnvr$LRR_median, mu = c(mu1, mu2, mu3),
                          mean.constr = c(mu1, NA, mu3), k = 3)
    
    paras_model$stage2$model <- list(mu = model3$mu,
                                     sigma = model3$sigma,
                                     lambda = model3$lambda)
    
    paras_all <- list(mus = model3$mu, 
                      sigmas = model3$sigma,
                      lambdas = model3$lambda)
    
  } else if (f1 == TRUE & f3 == FALSE) {
    
    mu3 <- mu2 + paras_LRR$LRR_mean$CN_3
    
    model3 <- normalmixEM(x = dt_cnvr$LRR_median, mu = c(mu1, mu2, mu3),
                          mean.constr = c(mu1, NA, mu3), k = 3)
    
    paras_model$stage2$model <- list(mu = model3$mu,
                                     sigma = model3$sigma,
                                     lambda = model3$lambda)
    
    paras_all <- list(mus = model3$mu, 
                      sigmas = model3$sigma,
                      lambdas = model3$lambda)
    
  } else {
    
    model3 <- normalmixEM(x = dt_cnvr$LRR_median, mu = c(mu1, mu2, mu3),
                          mean.constr = c(mu1, NA, mu3), k = 3)
    
    paras_model$stage2$model <- list(mu = model3$mu,
                                     sigma = model3$sigma,
                                     lambda = model3$lambda)
    
    paras_all <- list(mus = model3$mu, 
                      sigmas = model3$sigma,
                      lambdas = model3$lambda)
    
  }
  
  res <- list(paras_all = paras_all,
              paras_model = paras_model)
  return(res)  # return
  
}

model3_zz <- function(dt_cnvr, paras_LRR) {
  
  # init paras_model
  paras_model <- list()
  
  paras_model$stage1 <- list()
  paras_model$stage2 <- list()
  
  dt_cnvr <- dt_cnvr[order(dt_cnvr$CN), ]
  n <- nrow(dt_cnvr)  # number of samples
  numsnp <- unique(dt_cnvr$numSNP) # numsnp
  
  cn_factor <- factor(dt_cnvr$CN, levels = c(1, 2, 3))
  v <- as.vector(table(cn_factor))
  lambdas1 <- v/n ## 
  
  n1s <- v
  
  dt1 <- subset(dt_cnvr, CN == 1)
  dt2 <- subset(dt_cnvr, CN == 2)
  dt3 <- subset(dt_cnvr, CN == 3)
  
  mu3 <- mean(dt3$LRR_median)
  sigma3 <- sd(dt3$LRR_median)
  
  mu2 <- NA
  sigma2 <- NA
  
  flag1 <- TRUE  # 
  if (nrow(dt2) >= 200) { # 
    paras1 <- mode_zz(dt_cnvr = dt2)
    mu2 <- paras1$mu
    sigma2 <- paras1$sigma
  } else {
    paras1 <- mode_zz(dt_cnvr = dt_cnvr)
    mu2 <- paras1$mu
    
    mu3 <- mu2 + paras_LRR$LRR_mean$CN_3
    flag1 <- FALSE # changed flag1
  }
  
  # init paras
  mu1 <- mu2 + paras_LRR$LRR_mean$CN_1
  sigma1 <- paras_LRR$LRR_sd$CN_1
  
  paras_model$stage1$init <- list(mu = c(mu1, mu2, mu3),
                                  sigma = c(sigma1, sigma2, sigma3),
                                  lambda = lambdas1)

  model1 <- NULL
  flag2 <- TRUE
  if (flag1) {
    model1 <- normalmixEM(x = dt_cnvr$LRR_median, mu = c(mu1, mu2, mu3),
                          sigma = c(sigma1, sigma2, sigma3),
                          k = 3)
  } else {
    
    while(flag2) {
      
      model1 <- normalmixEM(x = dt_cnvr$LRR_median,
                            mu= c(mu1, mu2, mu3),
                            k = 3)
      mu1 <- model1$mu[1]
      mu2 <- model1$mu[2]
      mu3 <- model1$mu[3]
      
      d21 <- mu2 - mu1
      d32 <- mu3 - mu2
      
      f21 <- d21 >= 0.5*(-paras_LRR$LRR_mean$CN_1)
      f32 <- d32 >= 0.8*(paras_LRR$LRR_mean$CN_3)
      if (f21 & f32) {
        flag2 = FALSE
      }
      
    }
  }
  
  paras_model$stage1$model <- list(mu = model1$mu,
                                   sigma = model1$sigma,
                                   lambda = model1$lambda)
  
  if (flag2) {
    cns_pred <- apply(model1$posterior, MARGIN = 1, which.max)
    probs    <- model1$posterior
    # CN = 1
    probs1   <- as.vector(probs[, 1])
    idxs1    <- which(probs1 >= 0.9 & cns_pred == 1)
    n1 <- length(idxs1)
    n1
    
    # CN = 3
    probs3 <- as.vector(probs[, 3])
    idxs3 <- which(probs3 >= 0.9 & cns_pred == 3)
    n3 <- length(idxs3)
    n3
    
    idxs <- union(idxs1, idxs3)
    idxs2 <- setdiff(1:nrow(dt_cnvr), idxs)
    dt22 <- dt_cnvr[idxs2, ]
    n2 <- nrow(dt22)
    
    paras2 <- mode_zz(dt_cnvr = dt22)
    mu2 <- paras2$mu
    sigma2 <- paras2$sigma
    
    # parameters
    mu1 <- NA
    sigma1 <- NA
    if (n1 >= 10) {
      mu1 <- mean(dt_cnvr$LRR_median[idxs1])
      sigma1 <- sd(dt_cnvr$LRR_median[idxs1])
    } else {
      mu1 <- mu2 + paras_LRR$LRR_mean$CN_1
      sigma1 <- paras_LRR$LRR_sd$CN_1/sqrt(numsnp)
    }
    
    mu3 <- NA
    sigma3 <- NA
    if (n3 >= 10) {
      mu3 <- mean(dt_cnvr$LRR_median[idxs3])
      sigma3 <- sd(dt_cnvr$LRR_median[idxs3])
    } else {
      mu3 <- mu2 + paras_LRR$LRR_mean$CN_3
      sigma3 <- paras_LRR$LRR_sd$CN_3/sqrt(numsnp)
    }
    
    ## add
    if (mu3 < mu2) {
      mu3 <- mu2 + paras_LRR$LRR_mean$CN_3
      sigma3 <- paras_LRR$LRR_sd$CN_3/sqrt(numsnp)
    }
    
    if (mu1 > mu2) {
      mu1 <- mu2 + paras_LRR$LRR_mean$CN_1
      sigma1 <- paras_LRR$LRR_sd$CN_1/sqrt(numsnp)
    }
    
    n2 <- n - n1 - n3
    lambdas2 <- c(n1, n2, n3)/n
    
    n2s <- c(n1, n2, n3) 
    flag1 <- ifelse(n2s[1]/n1s[1] >= 0.8, TRUE, FALSE)
    flag3 <- ifelse(n2s[3]/n1s[3] >= 0.8, TRUE, FALSE)
    # init parameters
    paras_model$stage2$init <- list(mu = c(mu1, mu2, mu3),
                                    sigma = c(sigma1, sigma2, sigma3),
                                    lambda = lambdas2)
    
    model2 <- list()
    if (n1 >= 10 & n3 >= 10) {
      model2 <- normalmixEM(x = dt_cnvr$LRR_median, mu = c(mu1, mu2, mu3),
                            sigma = c(sigma1, sigma2, sigma3), 
                            mean.constr = c(mu1, NA, mu3),
                            k = 3)
    } else if (n1 >= 10 & n3 < 10) {
      model2 <- normalmixEM(x = dt_cnvr$LRR_median, mu = c(mu1, mu2, mu3),
                            sigma = c(sigma1, sigma2, sigma3), 
                            mean.constr = c(mu1, NA, mu3), k = 3)
    } else if (n1 < 10 & n3 >= 10) {
      model2 <- normalmixEM(x = dt_cnvr$LRR_median, mu = c(mu1, mu2, mu3),
                            sigma = c(sigma1, sigma2, sigma3), 
                            mean.constr = c(mu1, NA, mu3), k = 3)
    } else if (n1 < 10 & n3 < 10) {
      mu2 <- median(dt22$LRR_median)
      sigma2 <- sd(dt22$LRR_median)
      
      mu1 <- mu2 + paras_LRR$LRR_mean$CN_1
      if(n1 == 0) {
        sigma1 <- Inf
      } else {
        sigma1 <- paras_LRR$LRR_sd$CN_1/sqrt(numsnp)
      }
      
      mu3 <- mu2 + paras_LRR$LRR_mean$CN_3
      if (n3 == 0) {
        sigma3 <- Inf
      } else {
        sigma3 <- paras_LRR$LRR_sd$CN_3/sqrt(numsnp)
      }
      
      model2$mu <- c(mu1, mu2, mu3)
      model2$sigma <- c(sigma1, sigma2, sigma3)
      model2$lambda <- c(1, 1, 1)
    }
    
    # check parameters 
    cutoff1 <- 0.5
    cutoff3 <- 0.8
    mu1f <- model2$mu[1]
    mu2f <- model2$mu[2]
    mu3f <- model2$mu[3]
    
    f1 <- (abs(mu1f - mu2f) > cutoff1*abs(paras_LRR$LRR_mean$CN_1))
    f3 <- (abs(mu2f - mu3f) > cutoff3*abs(paras_LRR$LRR_mean$CN_3))
    
    f1 <- f1|flag1 ## add
    f3 <- f3|flag3
    
    if (f1 == FALSE & f3 == FALSE) {
      
      mu2 <- median(dt_cnvr$LRR_median)
      sd2 <- sd(dt_cnvr$LRR_median)
      
      mu3 <- mu2 + paras_LRR$LRR_mean$CN_3
      mu1 <- mu2 + paras_LRR$LRR_mean$CN_1
      
      model3 <- normalmixEM(x = dt_cnvr$LRR_median, mu = c(mu1, mu2, mu3), 
                            mean.constr = c(mu1, NA, mu3), k = 3)
      
      paras_model$stage2$model <- list(mu = model3$mu,
                                       sigma = model3$sigma,
                                       lambda = model3$lambda)
      
      paras_all <- list(mus = model3$mu, 
                        sigmas = model3$sigma,
                        lambdas = model3$lambda)
      
      
    } else if (f1 == FALSE & f3 == TRUE) {
      
      mu1 <- mu2 + paras_LRR$LRR_mean$CN_1
      
      model3 <- normalmixEM(x = dt_cnvr$LRR_median, mu = c(mu1, mu2, mu3),
                            mean.constr = c(mu1, NA, mu3), k = 3)
      
      paras_model$stage2$model <- list(mu = model3$mu,
                                       sigma = model3$sigma,
                                       lambda = model3$lambda)
      
      paras_all <- list(mus = model3$mu, 
                        sigmas = model3$sigma,
                        lambdas = model3$lambda)
      
    } else if (f1 == TRUE & f3 == FALSE) {
      
      mu3 <- mu2 + paras_LRR$LRR_mean$CN_3
      
      model3 <- normalmixEM(x = dt_cnvr$LRR_median, mu = c(mu1, mu2, mu3),
                            mean.constr = c(mu1, NA, mu3), k = 3)
      
      paras_model$stage2$model <- list(mu = model3$mu,
                                       sigma = model3$sigma,
                                       lambda = model3$lambda)
      
      paras_all <- list(mus = model3$mu, 
                        sigmas = model3$sigma,
                        lambdas = model3$lambda)
      
    } else {
      
      model3 <- normalmixEM(x = dt_cnvr$LRR_median, mu = c(mu1, mu2, mu3),
                            mean.constr = c(mu1, NA, mu3), k = 3)
      
      paras_model$stage2$model <- list(mu = model3$mu,
                                       sigma = model3$sigma,
                                       lambda = model3$lambda)
      
      paras_all <- list(mus = model3$mu, 
                        sigmas = model3$sigma,
                        lambdas = model3$lambda)
      
    }
    
    paras_model$stage2$model <- list(mu = model1$mu,
                                     sigma = model1$sigma,
                                     lambda = model1$lambda)
    res <- list(paras_all = paras_all,
                paras_model = paras_model)
    return(res)  # return
  } else {
    
    paras_model$stage2$init <- list(mu = model1$mu,
                                    sigma = model1$sigma,
                                    lambda = model1$lambda)
    paras_model$stage2$model <- list(mu = model1$mu,
                                     sigma = model1$sigma,
                                     lambda = model1$lambda)
    paras_all <- list(mus = model1$mu,
                      sigmas = model1$sigma,
                      lambdas = model1$lambda)
    
    res <- list(paras_all = paras_all,
                paras_model = paras_model)
    return(res)  # return
  }
  
}

model4_zz <- function(dt_cnvr, paras_LRR) {
  
  numsnp <- unique(dt_cnvr$numSNP)
  # init paras_model
  paras_model <- list()
  
  paras_model$stage1 <- list()
  paras_model$stage2 <- list()
  
  # only use one round gmm model
  dt_cnvr <- dt_cnvr[order(dt_cnvr$CN), ]
  n <- nrow(dt_cnvr)  # number of samples
  numsnp <- unique(dt_cnvr$numSNP) # numsnp
  
  cn_factor <- factor(dt_cnvr$CN, levels = c(1, 2, 3))
  ns <- as.vector(table(cn_factor))
  n1 <- ns[1]
  n2 <- ns[2]
  n3 <- ns[3]
  
  # start
  mu2 <- NA
  sigma2 <- NA
  
  dt1 <- subset(dt_cnvr, CN == 1)
  dt2 <- subset(dt_cnvr, CN == 2)
  dt3 <- subset(dt_cnvr, CN == 3)
  flag1 <- TRUE  # 
  if (nrow(dt2) >= 200) { # 
    paras1 <- mode_zz(dt_cnvr = dt2)
    mu2 <- paras1$mu
    sigma2 <- paras1$sigma
  } else {
    paras1 <- mode_zz(dt_cnvr = dt_cnvr)
    mu2 <- paras1$mu
    
    mu3 <- mu2 + paras_LRR$LRR_mean$CN_3
    flag1 <- FALSE # changed flag1
  }
  
  dt1 <- subset(dt_cnvr, CN == 1)
  mu1 <- mean(dt1$LRR_median)
  sigma1 <- sd(dt1$LRR_median)
  
  dt3 <- subset(dt_cnvr, CN == 3)
  mu3 <- mean(dt3$LRR_median)
  sigma3 <- sd(dt3$LRR_median)

  if (!flag1) {
    mu1 <- mu2 + paras_LRR$LRR_mean$CN_1
    sigma1 <- paras_LRR$LRR_sd$CN_1/sqrt(numsnp)
    
    mu3 <- mu2 + paras_LRR$LRR_mean$CN_3
    sigma3 <- paras_LRR$LRR_sd$CN_3/sqrt(numsnp)
  }
    
  
  
  # init paras
  paras_model$stage1$init <- list(mu = c(mu1, mu2, mu3),
                                  sigma = c(sigma1, sigma2, sigma3),
                                  lambda = c(n1, n2, n3)/n)
  
  
  model1 <- NULL
  flag2 <- TRUE
  if (flag1) {
    model1 <- normalmixEM(x = dt_cnvr$LRR_median, mu = c(mu1, mu2, mu3),
                          sigma = c(sigma1, sigma2, sigma3),
                          k = 3)  ## need to change
  } else {
    
    while(flag2) {
      
      model1 <- normalmixEM(x = dt_cnvr$LRR_median,
                            mu= c(mu1, mu2, mu3),
                            k = 3)
      mu1 <- model1$mu[1]
      mu2 <- model1$mu[2]
      mu3 <- model1$mu[3]
      
      d21 <- mu2 - mu1
      d32 <- mu3 - mu2
      
      f21 <- d21 >= 0.5*(-paras_LRR$LRR_mean$CN_1)
      f32 <- d32 >= 0.8*(paras_LRR$LRR_mean$CN_3)
      if (f21 & f32) {
        flag2 = FALSE
      }
      
    }
  }
  
  paras_model$stage1$model <- list(mu = model1$mu,
                                   sigma = model1$sigma,
                                   lambda = model1$lambda)
  
  if (flag2) {
    
    # number in CN = 1 and CN = 3
    cns_pred <- apply(model1$posterior, MARGIN = 1, which.max)
    probs <- model1$posterior
    
    probs1 <- as.vector(probs[, 1]) 
    idxs1  <- which(cns_pred == 1 & probs1 >= 0.99)
    n1 <- length(idxs1)
    
    probs3 <- as.vector(probs[, 3])
    idxs3  <- which(cns_pred == 3 & probs3 >= 0.99)
    n3 <- length(idxs3)
    
    idxs2 <- setdiff(1:nrow(dt_cnvr), union(idxs1, idxs3))
    dt22 <- dt_cnvr[idxs2, ]
    paras2 <- mode_zz(dt_cnvr = dt22)
    mu2 <- paras2$mu
    sigma2 <- paras2$sigma
    n2 <- nrow(dt22)
    
    mu1 <- 0
    sigma1 <- 0
    if (n1 >= 10) {
      dt21 <- dt_cnvr[idxs1, ]
      mu1 <- mean(dt21$LRR_median)
      sigma1 <- sd(dt21$LRR_median)
    } else {
      mu1 <- mu2 + paras_LRR$LRR_mean$CN_1
      sigma1 <- paras_LRR$LRR_sd$CN_1/sqrt(numsnp)
    }
    
    mu3 <- 0
    sigma3 <- 0
    if (n3 >= 10) {
      dt23 <- dt_cnvr[idxs3, ]
      mu3 <- mean(dt23$LRR_median)
      sigma3 <- sd(dt23$LRR_median)
    } else {
      mu3 <- mu2 + paras_LRR$LRR_mean$CN_3
      sigma3 <- paras_LRR$LRR_sd$CN_3/sqrt(numsnp)
    }
    
    ## add
    if (mu3 < mu2) {
      mu3 <- mu2 + paras_LRR$LRR_mean$CN_3
      sigma3 <- paras_LRR$LRR_sd$CN_3
    }
    
    n2 <- n - n1 - n3
    lambdas2 <- c(n1, n2, n3)/n
    # init parameters
    paras_model$stage2$init <- list(mu = c(mu1, mu2, mu3),
                                    sigma = c(sigma1, sigma2, sigma3),
                                    lambda = lambdas2)
    
    model2 <- list()
    if (n1 >= 10 & n3 >= 10) {
      model2 <- normalmixEM(x = dt_cnvr$LRR_median, mu = c(mu1, mu2, mu3),
                            sigma = c(sigma1, sigma2, sigma3), 
                            k = 3)
    } else if (n1 >= 10 & n3 < 10) {
      model2 <- normalmixEM(x = dt_cnvr$LRR_median, mu = c(mu1, mu2, mu3),
                            sigma = c(sigma1, sigma2, sigma3), 
                            mean.constr = c(NA, NA, mu3), k = 3)
    } else if (n1 < 10 & n3 >= 10) {
      model2 <- normalmixEM(x = dt_cnvr$LRR_median, mu = c(mu1, mu2, mu3),
                            sigma = c(sigma1, sigma2, sigma3), 
                            mean.constr = c(mu1, NA, mu3), k = 3)
    } else if (n1 < 10 & n3 < 10) {
      mu2 <- median(dt22$LRR_median)
      sigma2 <- sd(dt22$LRR_median)
      
      mu1 <- mu2 + paras_LRR$LRR_mean$CN_1
      mu3 <- mu2 + paras_LRR$LRR_mean$CN_3
      
      model5 <- normalmixEM(x = dt_cnvr$LRR_median, mu = c(mu1, mu2, mu3),
                            mean.constr = c(mu1, NA, mu3), k = 2)
      
      model2$mu <- model5$mu
      model2$sigma <- model5$sigma
      model2$lambda <- model5$lambda
    }
    
    # check parameters 
    cutoff1 <- 0.5
    cutoff3 <- 0.8
    mu1f <- model2$mu[1]
    mu2f <- model2$mu[2]
    mu3f <- model2$mu[3]
    
    f1 <- (abs(mu1f - mu2f) > cutoff1*abs(paras_LRR$LRR_mean$CN_1))
    f3 <- (abs(mu2f - mu3f) > cutoff3*abs(paras_LRR$LRR_mean$CN_3))
    
    if (f1 == FALSE & f3 == FALSE) {
      
      mu2 <- median(dt_cnvr$LRR_median)
      sd2 <- sd(dt_cnvr$LRR_median)
      
      mu3 <- mu2 + paras_LRR$LRR_mean$CN_3
      mu1 <- mu2 + paras_LRR$LRR_mean$CN_1
      
      model3 <- normalmixEM(x = dt_cnvr$LRR_median, mu = c(mu1, mu2, mu3), 
                            mean.constr = c(mu1, NA, mu3), k = 3)
      
      paras_model$stage2$model <- list(mu = model3$mu,
                                       sigma = model3$sigma,
                                       lambda = model3$lambda)
      
      paras_all <- list(mus = model3$mu, 
                        sigmas = model3$sigma,
                        lambdas = model3$lambda)
      
      
    } else if (f1 == FALSE & f3 == TRUE) {
      
      mu1 <- mu2 + paras_LRR$LRR_mean$CN_1
      
      model3 <- normalmixEM(x = dt_cnvr$LRR_median, mu = c(mu1, mu2, mu3),
                            mean.constr = c(mu1, NA, mu3), k = 3)
      
      paras_model$stage2$model <- list(mu = model3$mu,
                                       sigma = model3$sigma,
                                       lambda = model3$lambda)
      
      paras_all <- list(mus = model3$mu, 
                        sigmas = model3$sigma,
                        lambdas = model3$lambda)
      
    } else if (f1 == TRUE & f3 == FALSE) {
      
      mu3 <- mu2 + paras_LRR$LRR_mean$CN_3
      
      model3 <- normalmixEM(x = dt_cnvr$LRR_median, mu = c(mu1, mu2, mu3),
                            mean.constr = c(mu1, NA, mu3), k = 3)
      
      paras_model$stage2$model <- list(mu = model3$mu,
                                       sigma = model3$sigma,
                                       lambda = model3$lambda)
      
      paras_all <- list(mus = model3$mu, 
                        sigmas = model3$sigma,
                        lambdas = model3$lambda)
      
    } else {
      
      ## change mu1 mu2 mu3 value
      if (model2$mu[3] < model2$mu[2]) {
        model2$mu[3] = model2$mu[2] + paras_LRR$LRR_mean$CN_3
      } 
      if (model2$mu[1] > model2$mu[2]) {
        model2$mu[1] = model2$mu[2] - abs(paras_LRR$LRR_mean$CN_1)
      }
      
      paras_model$stage2$model <- list(mu = model2$mu,
                                       sigma = model2$sigma,
                                       lambda = model2$lambda)
      
      paras_all <- list(mus = model2$mu, 
                        sigmas = model2$sigma,
                        lambdas = model2$lambda)
      
    }
    
    res <- list(paras_all = paras_all,
                paras_model = paras_model)
    return(res)  # return
    
  } else {
    paras_model$stage2$init <- list(mu = model1$mu,
                                    sigma = model1$sigma,
                                    lambda = model1$lambda)
    paras_model$stage2$model <- list(mu = model1$mu,
                                     sigma = model1$sigma,
                                     lambda = model1$lambda)
    paras_all <- list(mus = model1$mu,
                      sigmas = model1$sigma,
                      lambdas = model1$lambda)
    
    res <- list(paras_all = paras_all,
                paras_model = paras_model)
    return(res)  # return
  }
  
}

# main function
train_model_zz <- function(dt_cnvr, paras_LRR) {
  
  cnvr_id <- unique(dt_cnvr$CNVR_ID)
  chr     <- unique(dt_cnvr$Chr)
  
  dt_cnvr <- dt_cnvr[order(dt_cnvr$CN), ]
  numsnp <- unique(dt_cnvr$numSNP)
  numsnp
  
  cns_factor <- factor(dt_cnvr$CN, levels = c(1, 2, 3)) 
  tbl <- table(cns_factor)
  tbl
  nums <- as.vector(tbl)
  nums
  n1 <- nums[1]
  n3 <- nums[3]
  
  # add permute data
  dt2 <- subset(dt_cnvr, CN == 2)
  mu2 <- NULL
  
  flag <- TRUE
  if (nrow(dt2) <= 100) { # search the mode for all data
    m2 <- mode_zz(dt_cnvr = dt_cnvr)
    mu2 <- m2$mu
    
    flag = FALSE
  } else {
    m2  <- mode_zz(dt_cnvr = dt2)
    mu2 <- m2$mu
  }
  
  mu1 <- mu2 + paras_LRR$LRR_mean$CN_1
  sigma1 <- paras_LRR$LRR_sd$CN_1/sqrt(numsnp)
  
  mu3 <- mu2 + paras_LRR$LRR_mean$CN_3
  sigma3 <- paras_LRR$LRR_sd$CN_3/sqrt(numsnp)
  
  num_permute <- 60 # 20
  
  LRR_median_cn1_permute <- rnorm(num_permute, mean = mu1, sd = sigma1)
  LRR_median_cn3_permute <- rnorm(num_permute, mean = mu3, sd = sigma3)
  
  dt_cn1_permute <- data.frame(Sample_ID = paste0("Sample_CN1_", 1:num_permute),
                               CNVR_ID = rep(cnvr_id, num_permute),
                               LRR_median = LRR_median_cn1_permute,
                               Chr = chr, alg = rep("permute", num_permute),
                               CN = 1, numSNP = numsnp, stringsAsFactors = FALSE)
  
  dt_cn3_permute <- data.frame(Sample_ID = paste0("Sample_CN3_", 1:num_permute),
                               CNVR_ID = rep(cnvr_id, num_permute),
                               LRR_median = LRR_median_cn3_permute,
                               Chr = chr, alg = rep("permute", num_permute),
                               CN = 3, numSNP = numsnp, stringsAsFactors = FALSE)
  
  # cutoff of minimum number of samples in each CN
  if (n1 >= 30 & n3 >= 30) {
    dt_cnvr <- dt_cnvr
    paras <- model1_zz(dt_cnvr = dt_cnvr, paras_LRR = paras_LRR)
    return(paras)
  } else if (n1 >= 30 & n3 < 30) {
    dt_cnvr <- rbind(dt_cnvr, dt_cn3_permute)
    paras <- model2_zz(dt_cnvr = dt_cnvr, paras_LRR = paras_LRR)
    return(paras)
  } else if (n1 < 30 & n3 >= 30) {
    
    if (flag == FALSE) {
      dt_cnvr <- rbind(dt_cnvr, dt_cn1_permute, dt_cn3_permute)
    } else {
      dt_cnvr <- rbind(dt_cnvr, dt_cn1_permute)
    }
    
    paras <- model3_zz(dt_cnvr = dt_cnvr, paras_LRR = paras_LRR)
    return(paras)
  } else if (n1 < 30 & n3 < 30) {
    dt_cnvr <- rbind(dt_cnvr, dt_cn1_permute, dt_cn3_permute)
    paras <- model4_zz(dt_cnvr = dt_cnvr, paras_LRR = paras_LRR)
    return(paras)
  }
  
}


