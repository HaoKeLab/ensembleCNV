

# process cnvrs LRR from snps to sample_based
process_cnvr_LRR <- function(dt_cnvrs, samples_LRR) {
  
  samples <- unique(dt_cnvrs$Sample_ID)
  # subset samples_LRR
  samples_LRR1 <- subset(samples_LRR, Sample_ID %in% samples)
  median_samples_LRR_SD <- median(samples_LRR1$LRR_SD, na.rm = TRUE)
  # test 
  stopifnot(nrow(samples_LRR1) == length(samples))
  
  res <- data.frame(Sample_ID = samples, CNVR_ID = unique(dt_cnvrs$CNVR_ID),
                    LRR_median = 0, Chr = unique(dt_cnvrs$Chr), alg = "other",
                    CN = 2, numSNP = unique(dt_cnvrs$numSNP), stringsAsFactors = FALSE)
  
  for (i in 1:length(samples)) {
    
    sample1 <- samples[i]
    
    idx1 <- which(samples_LRR1$Sample_ID == sample1)
    sample1_LRR_SD <- samples_LRR1$LRR_SD[idx1]  ##
    
    dt1 <- subset(dt_cnvrs, Sample_ID == sample1)
    
    LRR_median1 <- median(dt1$LRR, na.rm = TRUE)
    CN1 <- unique(dt1$CN)
    alg1 <- unique(dt1$alg)
    
    res$LRR_median[i] <- (LRR_median1/sample1_LRR_SD)*median_samples_LRR_SD ## transform
    # res$LRR_median[i] <- LRR_median1
    res$CN[i] <- CN1
    res$alg[i] <- alg1
  }
  
  res
}



# calculate LRR gatk whole with pi
calculate_LRR_gatk_whole <- function(dt_cnvr, mu1, sigma1, lambda1, cn_type) {
  if(cn_type == 2) { # for all CN = 2 type
    
    dt_cnvr$LRR2 <- sapply(1:nrow(dt_cnvr), FUN = function(k) {
      LRR1 <- dt_cnvr$LRR_median[k]
      prop1 <- lambda1*dnorm(x = LRR1, mean = mu1, sd = sigma1)
      prop1
    })
    
  } else if(cn_type == 1) {
    
    dt_cnvr$LRR1 <- sapply(1:nrow(dt_cnvr), FUN = function(k) {
      LRR1 <- dt_cnvr$LRR_median[k]
      prop1 <- lambda1*dnorm(x = LRR1, mean = mu1, sd = sigma1)
      prop1
    })
    
  } else if(cn_type == 3) {
    
    dt_cnvr$LRR3 <- sapply(1:nrow(dt_cnvr), FUN = function(k) {
      LRR1 <- dt_cnvr$LRR_median[k]
      prop1 <- lambda1*dnorm(x = LRR1, mean = mu1, sd = sigma1)
      prop1
    })
    
  } else if(cn_type == 0) {
    
    dt_cnvr$LRR0 <- sapply(1:nrow(dt_cnvr), FUN = function(k) {
      LRR1 <- dt_cnvr$LRR_median[k]
      prop1 <- lambda1*dnorm(x = LRR1, mean = mu1, sd = sigma1)
      prop1
    })
    
  }
  
  return(dt_cnvr)
}

# output LRR calcualte gatk result
output_LRR_gatk <- function(dt_cnvr, model) {
  
  dt_LRR0 <- calculate_LRR_gatk_whole(dt_cnvr = dt_cnvr, 
                                      mu1 = model$mu[1], 
                                      sigma1 = model$sigma[1], 
                                      lambda1 = model$lambda[1], cn_type = 0)
  dt_LRR1 <- calculate_LRR_gatk_whole(dt_cnvr = dt_cnvr, 
                                      mu1 = model$mu[2], 
                                      sigma1 = model$sigma[2], 
                                      lambda1 = model$lambda[2], cn_type = 1)
  dt_LRR2 <- calculate_LRR_gatk_whole(dt_cnvr = dt_cnvr, 
                                      mu1 = model$mu[3], 
                                      sigma1 = model$sigma[3], 
                                      lambda1 = model$lambda[3], cn_type = 2)
  dt_LRR3 <- calculate_LRR_gatk_whole(dt_cnvr = dt_cnvr, 
                                      mu1 = model$mu[4], 
                                      sigma1 = model$sigma[4], 
                                      lambda1 = model$lambda[4], cn_type = 3)
  dt_LRR01 <- merge(dt_LRR0, dt_LRR1)
  dt_LRR012 <- merge(dt_LRR01, dt_LRR2)
  dt_LRR0123 <- merge(dt_LRR012, dt_LRR3)  ## all p(LRR_median | CN = cn_type)
  
  return(dt_LRR0123)
}
