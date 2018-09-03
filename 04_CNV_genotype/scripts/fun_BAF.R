
# BAF emission probability (defined in PennCNV paper)
eBAF <- function (b, z, pB) {
  pib  <- 0.01
  
  mu0  <- 0.00
  mu14 <- 0.25 
  mu13 <- 1.0/3.0
  mu12 <- 0.5
  mu23 <- 2.0/3.0
  mu34 <- 0.75
  mu1  <- 1.00
  
  sd0  <- 0.016372
  sd14 <- 0.042099
  sd13 <- 0.045126
  sd12 <- 0.034982
  sd23 <- 0.045126
  sd34 <- 0.042099
  sd1  <- 0.016372
  
  M0   <- 0.5
  M1   <- 0.5
  
  sd5 <- 0.304243 ## for calculate CN = 0
  ## z=1, CN = 0, two copy deletion state
  if (z == 1) {
    e <- dnorm(b, mean = mu12, sd = sd5)
  }
  
  ## z=2, CN=1, one copy deletion state
  if (z==2) {
    e <- pib +
      (1 - pib) * (1-pB) * ( I(b==0)*M0 + I(b>0 & b<1)*(1-M0)*dnorm(b, mu0, sd0)/(1-pnorm(0,mu0,sd0)) ) +
      (1 - pib) * pB     * ( I(b==1)*M1 + I(b>0 & b<1)*(1-M1)*dnorm(b, mu1, sd1)/pnorm(1,mu1,sd1) )
  }
  
  ## z=3, CN=2, normal copy number state
  if (z==3) {
    e <- pib +
      (1 - pib) * 2*pB*(1-pB) * dnorm(b, mu12, sd12) +
      (1 - pib) * (1-pB)^2    * ( I(b==0)*M0 + I(b>0 & b<1)*(1-M0)*dnorm(b, mu0, sd0)/(1-pnorm(0,mu0,sd0)) ) +
      (1 - pib) * pB^2        * ( I(b==1)*M1 + I(b>0 & b<1)*(1-M1)*dnorm(b, mu1, sd1)/pnorm(1,mu1,sd1) )
  }
  
  ## z=4, CN=2, CN-LOH state
  if (z==4) {
    e <- pib +
      (1 - pib) * (1-pB) * ( I(b==0)*M0 + I(b>0 & b<1)*(1-M0)*dnorm(b, mu0, sd0)/(1-pnorm(0,mu0,sd0)) ) +
      (1 - pib) * pB     * ( I(b==1)*M1 + I(b>0 & b<1)*(1-M1)*dnorm(b, mu1, sd1)/pnorm(1,mu1,sd1) )
  }
  
  ## z=5, CN=3, one copy duplication state
  if (z==5) {
    e <- pib +
      (1 - pib) * 3*pB*(1-pB)^2 * dnorm(b, mu13, sd13) +
      (1 - pib) * 3*pB^2*(1-pB) * dnorm(b, mu23, sd23) +
      (1 - pib) * (1-pB)^3      * ( I(b==0)*M0 + I(b>0 & b<1)*(1-M0)*dnorm(b, mu0, sd0)/(1-pnorm(0,mu0,sd0)) ) +
      (1 - pib) * pB^3          * ( I(b==1)*M1 + I(b>0 & b<1)*(1-M1)*dnorm(b, mu1, sd1)/pnorm(1,mu1,sd1) )
  }
  
  return(e)
}

# BAF for gatk
baf_gatk_whole <- function(b, pB1, CN) {
  
  if (CN == 2) {
    return(eBAF(b = b, z = 3, pB = pB1))
  } else if (CN == 1) {
    return(eBAF(b = b, z = 2, pB = pB1))
  } else if (CN ==3) {
    return(eBAF(b = b, z = 5, pB = pB1))
  } else if (CN == 0) {
    return(eBAF(b = b, z = 1, pB = pB1))
  }
  
}

# calculate_BAF_gatk_whole CN = 0, 1, 2, 3
calculate_BAF_gatk_whole <- function(dt_cnvrs) {
  dt_cnvrs <- arrange(dt_cnvrs, Sample_ID, Name)
  samples <- unique(dt_cnvrs$Sample_ID)
  snps <- unique(dt_cnvrs$Name)
  snps <- dt_cnvrs$Name[1:length(snps)]  ## snps
  pfbs <- dt_cnvrs$PFB[1:length(snps)]  ## PFB
  
  m0 <- matrix(data = NA, nrow = length(samples), ncol = length(snps))
  m1 <- matrix(data = NA, nrow = length(samples), ncol = length(snps))
  m2 <- matrix(data = NA, nrow = length(samples), ncol = length(snps))
  m3 <- matrix(data = NA, nrow = length(samples), ncol = length(snps))
  
  for (i in 1:length(snps)) {
    
    snp1 <- snps[i]
    pfb1 <- pfbs[i]
    
    samples_snp <- subset(dt_cnvrs, Name == snp1)
    samples_new <- samples_snp$Sample_ID
    cns_new     <- samples_snp$CN  ## CN
    
    bafs_ep_0 <- sapply(samples_snp$BAF, FUN = function(x) baf_gatk_whole(b = x, pB1 = pfb1, CN = 0))
    bafs_ep_1 <- sapply(samples_snp$BAF, FUN = function(x) baf_gatk_whole(b = x, pB1 = pfb1, CN = 1))
    bafs_ep_2 <- sapply(samples_snp$BAF, FUN = function(x) baf_gatk_whole(b = x, pB1 = pfb1, CN = 2))
    bafs_ep_3 <- sapply(samples_snp$BAF, FUN = function(x) baf_gatk_whole(b = x, pB1 = pfb1, CN = 3))
    
    ## detect NaN values in bafs_ep_1/2/3
    idxs_na_0 <- which(is.na(bafs_ep_0))
    bafs_ep_0[idxs_na_0] <- median(bafs_ep_0, na.rm = TRUE)
    idxs_na_1 <- which(is.na(bafs_ep_1))
    bafs_ep_1[idxs_na_1] <- median(bafs_ep_1, na.rm = TRUE)  ## add median values
    idxs_na_2 <- which(is.na(bafs_ep_2))
    bafs_ep_2[idxs_na_2] <- median(bafs_ep_2, na.rm = TRUE)  ## add median values
    idxs_na_3 <- which(is.na(bafs_ep_3))
    bafs_ep_3[idxs_na_3] <- median(bafs_ep_3, na.rm = TRUE)  ## add median values
    
    m0[, i] <- bafs_ep_0
    m1[, i] <- bafs_ep_1
    m2[, i] <- bafs_ep_2
    m3[, i] <- bafs_ep_3
  }
  
  baf_eps_0 <- apply(m0, MARGIN = 1, prod)
  baf_eps_1 <- apply(m1, MARGIN = 1, prod) ### add na.rm
  baf_eps_2 <- apply(m2, MARGIN = 1, prod)
  baf_eps_3 <- apply(m3, MARGIN = 1, prod)
  
  dt_BAF <- data.frame(Sample_ID = samples_new, 
                       CN = cns_new, stringsAsFactors = FALSE)
  dt_BAF$BAF0 <- baf_eps_0
  dt_BAF$BAF1 <- baf_eps_1
  dt_BAF$BAF2 <- baf_eps_2
  dt_BAF$BAF3 <- baf_eps_3
  
  dt_BAF
}



