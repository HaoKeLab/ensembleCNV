

library(gridExtra)

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
    geom_histogram(bins = 50, fill = NA, color = "black") + 
    geom_line(data = dt123, aes(x, d, col = CN), lwd = 1.5) + 
    theme_bw(base_size = 10) + 
    labs(title = title,
         subtitle = paste("mu1:", round(mu1, 2), "mu2:", round(mu2, 2), "mu3:", round(mu3, 2), "\n",
                          "sd1:", round(sigma1, 2), "sd2:", round(sigma2, 2), "sd3:", round(sigma3, 2)))
  p
}

plot_gmm_diagnosis <- function(dt_cnvr, paras_model) {
  
  paras_stage1 <- paras_model$stage1
  paras_stage1_init <- paras_stage1$init
  paras_stage1_model <- paras_stage1$model
  
  paras_stage2 <- paras_model$stage2
  paras_stage2_init <- paras_stage2$init
  paras_stage2_model <- paras_stage2$model
  
  # plot 
  p1 <- plot_model(paras = paras_stage1_init, dt_cnvr = dt_cnvr, title = "stage1 init")
  p2 <- plot_model(paras = paras_stage1_model, dt_cnvr = dt_cnvr, title = "stage1 model")
    
  p3 <- plot_model(paras = paras_stage2_init, dt_cnvr = dt_cnvr, title = "stage2 init")
  p4 <- plot_model(paras = paras_stage2_model, dt_cnvr = dt_cnvr, title = "stage2 model")
    
  ps <- gridExtra::grid.arrange(p1, p2, p3, p4, nrow = 2)
  return(ps) 
}