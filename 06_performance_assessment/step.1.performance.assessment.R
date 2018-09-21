#!/usr/bin/env Rscript

suppressMessages(library(optparse))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))

option_list <- list(
  make_option(c("-d", "--duplicates"), action = "store", default = NA, type = "character",
              help = "Path to duplicate pairs information."),
  make_option(c("-n", "--matrixCN"), action = "store", default = NA,type = "character",
              help = "Path to matrix of copy number (CN)"),
  make_option(c("-g", "--matrixGQ"), action = "store", default = NA,type = "character",
              help = "Path to matrix of genotyping quality (GQ) score."),
  make_option(c("-o", "--resultpath"), action = "store", default = NA,type = "character",
              help = "Path to directory for saving assessment results.")
)

opt <- parse_args(OptionParser(option_list = option_list))
pars <- c(opt$duplicates, opt$matrixCN, opt$matrixGQ, opt$resultpath)

if (any(is.na(pars))) {
  stop("All required parameters must be supplied. (--help for detail)")
}

file_duplicates <- opt$duplicates
file_matrixcn   <- opt$matrixCN
file_matrixgq   <- opt$matrixGQ
path_result     <- opt$resultpath

dup_pairs <- read.delim(file = file_duplicates, as.is = TRUE)
matrix_CN <- readRDS(file = file_matrixcn)
matrix_gq <- readRDS(file = file_matrixgq)

# functions ---------------------------------------------------------------

generate_results <- function(mat, dup_pairs, gq1) {
  
  # clean cnvr with all CN = 2 or missing (denoted as -9)
  mat <- mat
  n_cnvr_raw <- nrow(mat)
  n_sample <- ncol(mat)
  
  # filter CNVR
  freqs <- unlist(lapply(1:n_cnvr_raw, FUN = function(k) {
    v1 <- as.vector(mat[k, ])
    freq1 <- sum(v1 %in% c(0, 1, 3))
    freq1
  }))
  
  idxs_del <- which(freqs == 0)
  mat1 <- mat
  if (length(idxs_del) >= 1) {
    mat1 <- mat[-idxs_del, ]
  }
  n_cnvr <- nrow(mat1)
  cat("After cleaning CNVRs with no CNV calls,", n_cnvr, "CNVRs remains from", n_cnvr_raw, "CNVRs.\n")
  mat_clean <- mat1  ## after cleaning nocall CNVR_ID
  
  ## sample level callRate
  freq_sample = unlist(lapply(1:n_sample, FUN = function(k) {
    v1 = as.vector( mat_clean[, k] )
    sum(v1 %in% c(0, 1, 2, 3))
  }))
  callRate_sample = freq_sample/n_cnvr
  ## CNVR level callRate
  freq_cnvr = unlist(lapply(1:n_cnvr, FUN = function(k) {
    v1 = as.vector( mat_clean[k, ] )
    sum(v1 %in% c(0, 1, 2, 3))
  }))
  callRate_cnvr = freq_cnvr/n_sample
  
  dat_callRate_cnvr = data.frame(callRate_cnvr = callRate_cnvr, 
                                 cutoff_gq = gq1, 
                                 stringsAsFactors = FALSE)
  dat_callRate_sample = data.frame(callRate_sample = callRate_sample,
                                   cutoff_gq = gq1, 
                                   stringsAsFactors = FALSE)
  
  ## consistency rate
  consistency_rates <- c()
  for (i in 1:nrow(dup_pairs)) {
    
    sample1 <- dup_pairs$sample1.name[i]
    sample2 <- dup_pairs$sample2.name[i]
    
    cns1 <- as.vector(mat1[, sample1])
    cns2 <- as.vector(mat1[, sample2])   # copy number of sample2
    
    # filter nocall cnvr
    idxs <- union(which(cns1 == -9), which(cns2 == -9))
    if(length(idxs) >= 1) {
      cns1 <- cns1[-idxs]
      cns2 <- cns2[-idxs]
    }
    
    idxs_overlap <- which(cns1 != 2 & cns2 != 2 & cns1 == cns2)
    idxs_union   <- union(which(cns1 != 2), which(cns2 != 2))
    
    rate1 <- length(idxs_overlap)/length(idxs_union)
    consistency_rates  <- c(consistency_rates, rate1)
  }
  
  res_consistency <- data.frame(consistency_rate = consistency_rates,
                                sample1.name = dup_pairs$sample1.name,
                                sample2.name = dup_pairs$sample2.name,
                                n_cnvr = n_cnvr,
                                cutoff_gq = gq1,
                                stringsAsFactors = FALSE)
  
  res_ncnvr = data.frame(cutoff_gq = gq1, n_cnvr = n_cnvr, 
                         stringsAsFactors = FALSE)
  ## return list results
  return(list(
    res_consistency = res_consistency,
    res_callRate_cnvr = dat_callRate_cnvr,
    res_callRate_sample = dat_callRate_sample,
    res_ncnvr = res_ncnvr
  ))
}

summary_regenotype = function(mat_CN, mat_gq, cutoffs_gq, dup_pairs) {
  
  res_consistency = data.frame() ## output dat for plot
  res_callRate_cnvr = data.frame()
  res_callRate_sample = data.frame()
  res_ncnvr = data.frame()
  
  for (i in 1:length(cutoffs_gq)) {
    
    gq1  = cutoffs_gq[i]
    cat("gq_score:", gq1, "\n")
    idx1 = which(mat_gq < gq1)
    if (length(idx1) >= 1) {
      mat_CN[idx1] = -9
    }
    
    res_gq1 = generate_results(mat = mat_CN,
                               dup_pairs = dup_pairs,
                               gq1 = gq1)
    
    res_consistency = rbind(res_consistency, res_gq1$res_consistency)
    res_callRate_sample = rbind(res_callRate_sample, res_gq1$res_callRate_sample)
    res_callRate_cnvr = rbind(res_callRate_cnvr, res_gq1$res_callRate_cnvr)
    res_ncnvr = rbind(res_ncnvr, res_gq1$res_ncnvr)
    
  }
  
  res_consistency$cutoff_gq = factor(res_consistency$cutoff_gq, levels = cutoffs_gq)
  res_callRate_sample$cutoff_gq = factor(res_callRate_sample$cutoff_gq, levels = cutoffs_gq)
  res_callRate_cnvr$cutoff_gq = factor(res_callRate_cnvr$cutoff_gq, levels = cutoffs_gq)
  res_ncnvr$cutoff_gq = factor(res_ncnvr$cutoff_gq, levels = cutoffs_gq)
  
  ## return list results
  return(list(
    res_consistency = res_consistency,
    res_callRate_cnvr = res_callRate_cnvr,
    res_callRate_sample = res_callRate_sample,
    res_ncnvr = res_ncnvr
  ))
}

# main --------------------------------------------------------------------

cutoffs_gq <- c(0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80)

res <- summary_regenotype(mat_CN = matrix_CN,
                          mat_gq = matrix_gq,
                          cutoffs_gq = cutoffs_gq,
                          dup_pairs = dup_pairs)

res_consistency     <- res$res_consistency
res_ncnvr           <- res$res_ncnvr
res_callRate_sample <- res$res_callRate_sample
res_callRate_cnvr   <- res$res_callRate_cnvr

saveRDS(res, file = file.path(path_result, "performance_assessment.rds"))

# start plot --------------------------------------------------------------

p1 = ggplot(data = res_consistency, aes(cutoff_gq, consistency_rate)) +
  geom_boxplot() +
  theme_bw() + 
  geom_hline(yintercept = 0.9, lty = 2, lwd = 1, col = "grey60") +
  theme(panel.background = element_blank(),
        plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text = element_text(size = 15),
        strip.text = element_text(size = 15)) +
  xlab("GQ score threshold") +
  ylab("Concordance rate") +
  ggtitle("Concrodance rate")

p2 = ggplot(data = res_ncnvr, aes(cutoff_gq, n_cnvr)) +
  #geom_col() +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(panel.background = element_blank(),
        plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text = element_text(size = 15),
        strip.text = element_text(size = 15)) +
  xlab("GQ score threshold") +
  ylab("Number of CNVRs") +
  ggtitle("Number of CNVRs")

p3 = ggplot(data = res_callRate_sample, aes(cutoff_gq, callRate_sample)) +
  geom_boxplot() +
  theme_bw() +
  geom_hline(yintercept = 0.9, lty = 2, lwd = 1, col = "grey60") +
  theme(panel.background = element_blank(),
        plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text = element_text(size = 15),
        strip.text = element_text(size = 15)) +
  xlab("GQ score threshold") +
  ylab("Sample-wise call rate") +
  ggtitle("Sample-wise call rate")


p4 = ggplot(data = res_callRate_cnvr, aes(cutoff_gq, callRate_cnvr)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.9, lty = 2, lwd = 1, col = "grey60") +
  theme_bw() +
  theme(panel.background = element_blank(),
        plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text = element_text(size = 15),
        strip.text = element_text(size = 15)) +
  xlab("GQ score threshold") +
  ylab("CNVR-wise call rate") +
  ggtitle("CNVR-wise call rate")

png(filename = file.path(path_result, "performance_assessment.png"),
    width = 12, height = 12, units = "in", res = 512)

p = plot_grid(p1, p2, p3, p4, 
              nrow = 2, ncol = 2,
              labels = LETTERS[1:4],
              label_size = 22,
              vjust = 1.2, align = "hv")
print(p)

dev.off()


