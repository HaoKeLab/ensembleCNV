

# add dup pairs flag
add_dup_pairs_flag <- function(dt, dup_pairs) {
  
  # dup pair with dup_flag not equal to 0
  dt$dup_flag <- 0
  for (i in 1:nrow(dup_pairs)) {
    samples <- c(dup_pairs$sample1.name[i], dup_pairs$sample2.name[i])
    idxs <- which(dt$Sample_ID %in% samples)
    if (length(idxs) >= 1) {
      dt$dup_flag[idxs] <- i
    }
  }
  
  dt
}

# tranfrom dt_LRRBAF from gatk to raw method we used LRR12/BAF12 
add_LRRBAF_ratio <- function(dt) {
  
  # add log_ratio LRR12 and BAF12 and LRR32 and BAF32
  dt$LRR12 <- log(dt$LRR1/dt$LRR2)
  dt$BAF12 <- log(dt$BAF1/dt$BAF2)
  
  dt$LRR32 <- log(dt$LRR3/dt$LRR2)
  dt$BAF32 <- log(dt$BAF3/dt$BAF2)
  
  return(dt)
}

# plot model
plot_model <- function(paras, dt_cnvr, title) {
  
  mu1 <- paras$mu[1]
  sigma1 <- paras$sigma[1]
  lambda1 <- paras$lambda[1]
  
  mu2 <- paras$mu[2]
  sigma2 <- paras$sigma[2]
  lambda2 <- paras$lambda[2]
  
  mu3 <- paras$mu[3]
  sigma3 <- paras$sigma[3]
  lambda3 <- paras$lambda[3]
  
  x <- dt_cnvr$LRR_median
  range_x <- range(x)
  
  xs <- seq(range_x[1], range_x[2], length.out = 800)
  dt <- data.frame(x = xs, stringsAsFactors = F)
  
  dt1 <- data.frame(x = xs, d = lambda1*dnorm(xs, mean = mu1, sd = sigma1), CN = 1)
  dt3 <- data.frame(x = xs, d = lambda3*dnorm(xs, mean = mu3, sd = sigma3), CN = 3)
  dt2 <- data.frame(x = xs, d = lambda2*dnorm(xs, mean = mu2, sd = sigma2), CN = 2)
  dt123 <- rbind(dt1, dt2, dt3)
  dt123$CN <- as.factor(dt123$CN)
  
  p <- ggplot(data = dt_cnvr, aes(LRR_median, y = ..density..)) +
    geom_histogram(bins = 100, fill = NA, color = "black") + 
    geom_line(data = dt123, aes(x, d, col = CN), lwd = 1.5) + 
    theme_bw(base_size = 10) +
    labs(title = title,
         subtitle = paste("mu1:", round(mu1, 2), "mu2:", round(mu2, 2), "mu3:", round(mu3, 2), "\n",
                          "sd1:", round(sigma1, 2), "sd2:", round(sigma2, 2), "sd3:", round(sigma3, 2)))
  p
}


# plot steps
plot_steps <- function(dt_cnvr_train, dup_pairs, paras, dt_cnvr_raw, dt_LRRBAF) {
  
  dt_cnvr_train <- dt_cnvr_train[order(dt_cnvr_train$CN), ]
  dt_cnvr_train$idx <- 1:nrow(dt_cnvr_train)
  
  dt_cnvr_raw <- dt_cnvr_raw[order(dt_cnvr_raw$CN), ]
  dt_cnvr_raw$idx <- 1:nrow(dt_cnvr_raw)
  
  dt_dup     <- data.frame()
  dt_dup_raw <- data.frame()
  if (! is.null(dup_pairs) ) {
    # add flag for dup
    dt_cnvr_train <- add_dup_pairs_flag(dt = dt_cnvr_train, dup_pairs = dup_pairs)
    dt_dup <- subset(dt_cnvr_train, dup_flag != 0)
    
    dt_cnvr_raw <- add_dup_pairs_flag(dt = dt_cnvr_raw, dup_pairs = dup_pairs)
    dt_dup_raw <- subset(dt_cnvr_raw, dup_flag != 0)
  }
  
  numsnp <- unique(dt_cnvr_raw$numSNP)
  
  plot1 <- ggplot(data = dt_cnvr_raw, aes(idx, LRR_median, col = factor(CN))) + 
    geom_point() + 
    theme_bw(base_size = 10) +
    annotate("text", x = dt_dup_raw$idx, y = dt_dup_raw$LRR_median, label = dt_dup_raw$dup_flag) + 
    labs(title = paste("scatter plot of LRR_median with numsnp:", numsnp)) + 
    theme(legend.position = "top")
  
  # plot1
  
  # add gmm model paras
  plot2 <- plot_model(dt_cnvr = dt_cnvr_train, paras = paras, 
                      title = "fit model for LRR_median, only contain CN = 1/2/3")
  
  # plot2
  
  # --------------------------------------------------
  # add steps infromation 
  # filter CN != 0
  dt_LRRBAF <- subset(dt_LRRBAF, CN_gatk_pred != 0)
  
  dt_LRRBAF <- add_dup_pairs_flag(dt = dt_LRRBAF, dup_pairs = dup_pairs)
  dt_LRRBAF_new <- add_LRRBAF_ratio(dt = dt_LRRBAF)
  
  # step1
  dt1_gatk <- subset(dt_LRRBAF_new, CN_gatk_pred == 1)
  dt1_annotate <- subset(dt_LRRBAF_new, dup_flag != 0)
  plot_step1 <- ggplot() + 
    geom_point(data = dt_LRRBAF_new, aes(BAF12, LRR12), col = "gray") + 
    geom_vline(xintercept = 0, lty = 2, lwd = 1) + 
    geom_hline(yintercept = 0, lty = 2, lwd = 1) + 
    theme_bw(base_size = 10) +
    geom_point(data = dt1_gatk, aes(BAF12, LRR12), col = "red") + 
    annotate(geom = "text", x = dt1_annotate$BAF12, y = dt1_annotate$LRR12, label = dt1_annotate$dup_flag) + 
    ggtitle(label = "step 1 for CN = 1")
  
  # step2
  dt3 <- subset(dt_LRRBAF_new, LRR12 <= 0 | BAF12 <= 0)
  dt3_gatk <- subset(dt_LRRBAF_new, CN_gatk_pred == 3)
  dt3_annotate <- subset(dt3, dup_flag != 0)
  plot_step2 <- ggplot() + 
    geom_point(data = dt3, aes(BAF32, LRR32), col = "gray") + 
    geom_vline(xintercept = 0, lty = 2, lwd = 1) + 
    geom_hline(yintercept = 0, lty = 2, lwd = 1) + 
    theme_bw(base_size = 10) +
    geom_point(data = dt3_gatk, aes(BAF32, LRR32), col = "red") + 
    annotate(geom = "text", x = dt3_annotate$BAF32, y = dt3_annotate$LRR32, label = dt3_annotate$dup_flag) + 
    ggtitle(label = "step 2 for CN = 3")
  
  
  ps <- gridExtra::grid.arrange(plot1, plot_step1, plot2, plot_step2, nrow = 2)
  
  return(ps)
}


# plot model_final
plot_model_final <- function(paras, dt_cnvr, title) {
  
  mu1 <- paras$mu[2]
  sigma1 <- paras$sigma[2]
  lambda1 <- paras$lambda[2]
  
  mu2 <- paras$mu[3]
  sigma2 <- paras$sigma[3]
  lambda2 <- paras$lambda[3]
  
  mu3 <- paras$mu[4]
  sigma3 <- paras$sigma[4]
  lambda3 <- paras$lambda[4]
  
  
  # transfrom lambdas --------------
  lambdas <- lambda1 + lambda2 + lambda3
  lambda1 <- lambda1/lambdas
  lambda2 <- lambda2/lambdas
  lambda3 <- lambda3/lambdas
  
  x <- dt_cnvr$LRR_median
  range_x <- range(x)
  
  xs <- seq(range_x[1], range_x[2], length.out = 800)
  dt <- data.frame(x = xs, stringsAsFactors = F)
  
  dt1 <- data.frame(x = xs, d = lambda1*dnorm(xs, mean = mu1, sd = sigma1), CN = 1)
  dt2 <- data.frame(x = xs, d = lambda2*dnorm(xs, mean = mu2, sd = sigma2), CN = 2)
  dt3 <- data.frame(x = xs, d = lambda3*dnorm(xs, mean = mu3, sd = sigma3), CN = 3)
 
  dt123 <- rbind(dt1, dt2, dt3)
  dt123$CN <- as.factor(dt123$CN)
  
  p <- ggplot(data = dt_cnvr, aes(LRR_median, y = ..density..)) +
    geom_histogram(bins = 100, fill = NA, color = "black") + 
    geom_line(data = dt123, aes(x, d, col = CN), lwd = 1.5) + 
    theme_bw(base_size = 10) +
    labs(title = title,
         subtitle = paste("mu1:", round(mu1, 2), "mu2:", round(mu2, 2), "mu3:", round(mu3, 2), "\n",
                          "sd1:", round(sigma1, 2), "sd2:", round(sigma2, 2), "sd3:", round(sigma3, 2)))
  p
}
