#!/usr/bin/env Rscript --vanilla

args <- commandArgs(trailingOnly = TRUE)

path_input  <- args[1]
cohort_name <- args[2]
path_output <- args[3]

# data --------------------------------------------------------------------

compact_info <- function(m1, sample1, sample2, type) {
  
  cns_1 <- as.integer(m1[, sample1])
  cns_2 <- as.integer(m1[, sample2])
  
  ## del nocall CNVRs (-9)
  idx.nocall = union( which(cns_1 == -9), which(cns_2 == -9) )
  if ( length(idx.nocall) >= 1 ) {
    cns_1 = cns_1[-idx.nocall]
    cns_2 = cns_2[-idx.nocall]
  } 
  
  idxs_overlap <- which(cns_1 == cns_2 & cns_1 != 2 & cns_2 != 2)
  idxs_union <- union(which(cns_1 != 2), which(cns_2 != 2))
  
  dt1 <- data.frame(sample1 = sample1, sample2 = sample2,
                    prop_overlap = length(idxs_overlap)/length(idxs_union),
                    num_union = length(idxs_union), num_overlap = length(idxs_overlap),
                    cnv_method = type, stringsAsFactors = FALSE)
  dt1
}

prepare_data <- function(path) {
  
  dup_pairs = readRDS(file = file.path(path, "dup_samples.rds"))
  
  mat_iPattern <- readRDS(file = file.path(path, "matrix_iPattern.rds"))
  mat_PennCNV <- readRDS(file = file.path(path, "matrix_PennCNV.rds"))
  mat_QuantiSNP <- readRDS(file = file.path(path, "matrix_QuantiSNP.rds"))
  mat_ensembleCNV_intersect <- readRDS(file = file.path(path, "matrix_IPQ_intersect.rds"))
  mat_ensembleCNV_union <- readRDS(file = file.path(path, "matrix_IPQ_union.rds"))
  mat_ensembleCNV <- readRDS(file = file.path(path, "matrix_ensembleCNV.rds"))
  
  n.dups  <- nrow(dup_pairs)
  n.types <- 6
  types   <- c("iPattern", "PennCNV", "QuantiSNP", "intersect", "union", "ensembleCNV")
  
  res_compact <- data.frame(sample1 = rep(dup_pairs$sample1.name, n.types),
                            sample2 = rep(dup_pairs$sample2.name, n.types),
                            prop_overlap = rep(0, n.dups*n.types),
                            num_union = rep(0, n.dups*n.types),
                            num_overlap = rep(0, n.dups*n.types),
                            LRR_type = rep(types, each = n.dups), stringsAsFactors = FALSE)
  
  for(i in 1:nrow(dup_pairs)) {
    
    sample1 <- dup_pairs$sample1.name[i]
    sample2 <- dup_pairs$sample2.name[i]
    
    cat(i, "in", nrow(dup_pairs),"dup_pair:", sample1, sample2, "\n")
    # for compact_info
    res_compact[i, ] <- compact_info(m1 = mat_iPattern, sample1 = sample1, sample2 = sample2, type = "iPattern")
    res_compact[(i+n.dups), ] <- compact_info(m1 = mat_PennCNV, sample1 = sample1, sample2 = sample2, type = "PennCNV")
    res_compact[(i+2*n.dups), ] <- compact_info(m1 = mat_QuantiSNP, sample1 = sample1, sample2 = sample2, type = "QuantiSNP")
    res_compact[(i+3*n.dups), ] <- compact_info(m1 = mat_ensembleCNV_intersect, sample1 = sample1, sample2 = sample2, type = "intersect")
    res_compact[(i+4*n.dups), ] <- compact_info(m1 = mat_ensembleCNV_union, sample1 = sample1, sample2 = sample2, type = "union")
    res_compact[(i+5*n.dups), ] <- compact_info(m1 = mat_ensembleCNV, sample1 = sample1, sample2 = sample2, type = "ensembleCNV")
  }
  
  res_compact$prop_overlap <- as.numeric(res_compact$prop_overlap)
  res_compact <- na.omit(res_compact)
  
  res_compact$LRR_type <- factor(res_compact$LRR_type, levels = types)
  
  ## callRate sample level
  n_sample = ncol(mat_ensembleCNV)
  n_cnvr = nrow(mat_ensembleCNV)
  
  callRate_sample = unlist(lapply(1:n_sample, FUN = function(k) {
    v1 = as.vector(mat_ensembleCNV[, k])
    sum(v1 %in% c(0, 1, 2, 3))/n_cnvr
  }))
  
  callRate_cnvr = unlist(lapply(1:n_cnvr, FUN = function(k) {
    v1 = as.vector(mat_ensembleCNV[k, ])
    sum(v1 %in% c(0, 1, 2, 3))/n_sample
  }))
  
  dat_callRate_Sample = data.frame(callRate_Sample = callRate_sample, stringsAsFactors = FALSE)
  dat_callRate_CNVR = data.frame(callRate_CNVR = callRate_cnvr, stringsAsFactors = FALSE)
  
  
  return(list(
    res_compact = res_compact,
    dat_callRate_Sample = dat_callRate_Sample,
    dat_callRate_CNVR = dat_callRate_CNVR
  ))
}

# copy all following files to path_input:
# dup_samples.rds with columns: sample1.name, sample1.name
# matrix_iPattern.rds; matrix_PennCNV.rds; matrix_QuantiSNP.rds; 
# matrix_IPQ_intersect.rds; matrix_IPQ_union.rds; matrix_ensembleCNV.rds

data_plot <- prepare_data(path = path_input)
saveRDS(data_plot, file = file.path( path_output, "dups.consistency.rds"))

# plot --------------------------------------------------------------------

plots <- function(data, cohort) {
  
  require(ggplot2, quietly = TRUE)
  
  data$res_compact$LRR_type <- as.character(data$res_compact$LRR_type)
  data$res_compact$LRR_type[which(data$res_compact$LRR_type == "intersect")] <- "intersection"
  data$res_compact$LRR_type <- factor(data$res_compact$LRR_type, 
                                      levels = c("iPattern", "PennCNV", "QuantiSNP", "intersection", "union","ensembleCNV"))
  p1 <- ggplot(data = data$res_compact, 
               aes(LRR_type, prop_overlap, fill = LRR_type)) + 
    geom_boxplot() + 
    theme_bw(base_size = 9) + 
    theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
          axis.text.x = element_text(size = 12, angle = 15, vjust = 1, hjust = 1, face = "bold"),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 13, margin = margin(r = 5), face = "bold"),
          legend.position = "none") +
    ggtitle(paste0(cohort, ": concordance rate")) +
    xlab("") + 
    ylab("Concordance rate") + 
    scale_y_continuous(breaks = seq(0, 1, 0.1),
                       limits = c(0, 1))
  
  dat1 = ggplot_build(p1)$data[[1]]
  col1 = dat1$fill[6]
  
  p2 <- ggplot(data = data$dat_callRate_Sample, aes(callRate_Sample)) + 
    geom_histogram(bins = 50, fill = col1) +
    theme_bw(base_size = 9) +
    theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 13, margin = margin(r = 5), face = "bold"),
          axis.title.x = element_text(size = 13, margin = margin(t = -35), face = "bold")) + 
    ggtitle(paste0(cohort, ": sample-wise call rate")) +
    xlab("Call rate") +
    ylab("Frequency") +
    scale_x_continuous(limits = c(0.5, 1))
  
  p3 = ggplot(data = data$dat_callRate_CNVR, aes(callRate_CNVR, fill = col1)) + 
    geom_histogram(bins = 50,fill = col1) + 
    theme_bw(base_size = 9) + 
    theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 13, margin = margin(r = 5), face = "bold"),
          axis.title.x = element_text(size = 13, margin = margin(t = -35), face = "bold")) + 
    ggtitle(paste0(cohort, ": CNVR-wise call rate")) +
    xlab("Call rate") +
    ylab("Frequency")
  
  return(list(
    p1 = p1,
    p2 = p2,
    p3 = p3
  ))
}

require(ggplot2)
require(cowplot)

plots_all <- plots(data = data_plot, cohort = cohort_name)

png( filename = file.path(path_output, "dups.consistency.png"),
     width = 18, height = 6, units = "in", res = 512)

pg <- plot_grid(plots_all$p1, plots_all$p2, plots_all$p3,
          nrow = 1, ncol = 3, labels = LETTERS[1:6], label_size = 20,
          vjust = 1.2, align = "h")
print(pg)

dev.off()



