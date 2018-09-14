#!/usr/bin/env Rscript

suppressMessages({
  library(pheatmap)
  library(Rcpp)
  library(optparse)
})

option_list <- list(
  make_option(c("-c", "--chr"), action = "store", type = "integer", default = NA,
              help = "Select CNVRs on chromosome CHR for boundary refinement."),
  make_option(c("-p", "--datapath"), action = "store", type = "character", default = NA,
              help = "Path to the directory containing necessary input data."),
  make_option(c("-o", "--resultpath"), action = "store", type = "character", default = NA,
              help = "Path to the directory for saving results."),              
  make_option(c("-m", "--matrixpath"), action = "store", type = "character", default = NA, 
              help = "Path to chromosome-wise LRR and BAF matrices."),
  make_option(c("-s", "--rcppfile"), action = "store", type = "character", default = NA,
              help = "Path to refine.rcpp to be used in this R script."),
  make_option(c("-r", "--centromere"), action = "store", type = "character", default = NA,
              help = "Path to file with centromere position of each chromosome."),
  make_option(c("-n", "--plot"), action = "store_true", default = FALSE,
              help = "[optional] Whether to generate diagnosis plots.")
)

opt <- parse_args(OptionParser(option_list = option_list))
pars = c(opt$chr, opt$datapath, opt$resultpath, opt$matrixpath, opt$rcppfile, opt$centromere)

if ( any(is.na(pars)) ) {
  stop("All parameters must be supplied. (--help for detail)")
}

chr1        <- opt$chr
path_data   <- opt$datapath
path_matrix <- opt$matrixpath
path_result <- opt$resultpath ## cnvr_refine.txt located in this directory
file_rcpp   <- opt$rcppfile
flag_plot   <- opt$plot

# cnvrs refinement
dat_cnvrs_refine <- read.delim( file = file.path(path_result, "cnvr_refine.txt"), as.is = TRUE ) 
stopifnot( nrow(dat_cnvrs_refine) > 0 )

dat_cnvrs_refine_chr1 <- subset(dat_cnvrs_refine, chr == chr1)
stopifnot( nrow(dat_cnvrs_refine_chr1) > 0 )

# LRR matrix
file_LRR <- paste0("matrix_chr_", chr1, "_LRR.rds")
dt_matrix_LRR <- readRDS(file = file.path(path_matrix, "LRR", file_LRR))
mat_LRR <- as.matrix(dt_matrix_LRR)

samples_LRR <- rownames( mat_LRR )
snps_LRR <- colnames( mat_LRR )
n_snp <- ncol(mat_LRR)

# fill NA values in mat_LRR
snps_mean <- colMeans( mat_LRR, na.rm = T)
for ( i in 1:nrow(mat_LRR) ) {
  v1   <- mat_LRR[i, ]
  idx1 <- which(is.na(v1))
  if (length(idx1) >= 1) {
    mat_LRR[i, idx1] <- snps_mean[idx1]
  }
}

# centromere position (hg19)
centromere <- read.delim(file = opt$centromere, as.is = TRUE)
pos_centromere_chr1 <- centromere$position[centromere$chr == chr1]

# position
dat_pfb <- read.table(file = file.path(path_data, "SNP.pfb"), sep = "\t",
                      header = TRUE, as.is = TRUE, check.names = FALSE,
                      comment.char = "")
dat_pfb_chr1 <- subset(dat_pfb, Chr == chr1)
snp_chr1 <- dat_pfb_chr1

snps_chr1_p <- subset(dat_pfb_chr1, Position <= pos_centromere_chr1)
snps_chr1_q <- subset(dat_pfb_chr1, Position > pos_centromere_chr1)

snps_chr1_p <- snps_chr1_p[order(snps_chr1_p$Position, snps_chr1_p$Name), ]
snps_chr1_q <- snps_chr1_q[order(snps_chr1_q$Position, snps_chr1_q$Name), ]

mat_LRR_p <- mat_LRR[, snps_chr1_p$Name]
mat_LRR_q <- mat_LRR[, snps_chr1_q$Name]

snps_LRR_p <- colnames( mat_LRR_p )
snps_LRR_q <- colnames( mat_LRR_q )

n_snps_p <- length( snps_LRR_p )
n_snps_q <- length( snps_LRR_q )

# rcpp create CNVR --------------------------------------------------------
sourceCpp(file = file_rcpp)

# main --------------------------------------------------------------------
path_refine <- file.path(path_result, "res_refine")
if (!dir.exists(paths = path_refine)) dir.create(path = path_refine, showWarnings = F, recursive = T)

folder_chr1 <- paste0("chr", chr1)
path_png <- file.path(path_refine, folder_chr1, "png")
path_res <- file.path(path_refine, folder_chr1, "data")

dir.create(path = path_png, showWarnings = F, recursive = T)
dir.create(path = path_res, showWarnings = F, recursive = T)

res <- data.frame()
n_cnvrs_chr1 <- nrow(dat_cnvrs_refine_chr1)
cat("number of CNVRs to be boundary-refined:", n_cnvrs_chr1, "\n")

for ( i in 1:n_cnvrs_chr1 ) {
  
  cnvr1     <- dat_cnvrs_refine_chr1$CNVR_ID[i]
  chr1      <- dat_cnvrs_refine_chr1$chr[i]
  snp_start <- dat_cnvrs_refine_chr1$start_snp[i]
  snp_end   <- dat_cnvrs_refine_chr1$end_snp[i]
  freq1     <- dat_cnvrs_refine_chr1$Freq[i]
  
  cat(i, "in", n_cnvrs_chr1, "CNVR_ID:", cnvr1, "\n")
  # p or q arm
  strs <- unlist( strsplit(cnvr1, split = "_", fixed = TRUE))
  arm_cnvr1 <- strs[ length(strs) ]
  
  mat_LRR_cnvr1 <- NULL
  snp_LRR_cnvr1 <- NULL
  n_snps_cnvr1  <- NULL
  if ( arm_cnvr1 == "p" ) {
    mat_LRR_cnvr1 <- mat_LRR_p
    snp_LRR_cnvr1 <- snps_LRR_p
    n_snps_cnvr1  <- n_snps_p
  } else if ( arm_cnvr1 == "q" ) {
    mat_LRR_cnvr1 <- mat_LRR_q
    snp_LRR_cnvr1 <- snps_LRR_q
    n_snps_cnvr1  <- n_snps_q
  } else {
    stop(paste(cnvr1, 'must be in the p or q arm.'))
  }
  # raw CNVR information
  idx_start = which(snp_LRR_cnvr1 == snp_start)
  idx_end = which(snp_LRR_cnvr1 == snp_end)
  stopifnot( idx_end > idx_start )
  n_snp_chr1 = idx_end - idx_start + 1
  
  # check number of SNPs in the cnvr1
  # set cutoff n.snps <= 100
  if ( n_snp_chr1 >= 100 ) {
    
    snp.pos.start = dat_pfb_chr1$Position[which(dat_pfb_chr1$Name == snp_start)]
    snp.pos.end = dat_pfb_chr1$Position[which(dat_pfb_chr1$Name == snp_end)] 
    
    # plot pheatmap
    idx11 = which(snp_LRR_cnvr1 == snp_start)
    idx22 = which(snp_LRR_cnvr1 == snp_end)
    
    idx11.new = ifelse((idx11 - 10) <= 0, 1, idx11 - 10)
    idx22.new = ifelse((idx22+10) >= n_snps_cnvr1, n_snps_cnvr1, idx22+10)
    
    mat_cnvr1_raw = mat_LRR_cnvr1[, idx11.new:idx22.new]
    
    len1 = ifelse((idx11-10) <= 0, idx11-1, 10)
    len2 = ifelse((idx22+10) >= n_snps_cnvr1, n_snps_cnvr1 - idx22, 10)
    
    snps_flag_raw = c(rep("out", len1), rep("in", n_snp_chr1), rep("out", len2))
    
    annotate_col1 = data.frame(
      group_raw = snps_flag_raw
    )
    
    mcorr_raw = cor(mat_cnvr1_raw)
    rownames(annotate_col1) = colnames(mat_cnvr1_raw)
    
    if ( flag_plot ) {
      filename_png = paste0(cnvr1, "_boundary_refinement.png")
      png(filename = file.path(path_png, filename_png), width = 12, height = 12, res = 512, units = "in")
      par(mar = c(4, 4, 4, 4))
      pheatmap(mcorr_raw,
               cluster_rows = FALSE,
               cluster_cols = FALSE,
               annotation_col = annotate_col1,
               show_rownames = FALSE,
               show_colnames = FALSE,
               main = paste(cnvr1, "raw & not do refine"))
      dev.off()
    }
    
    res1 = data.frame(CNVR_ID = cnvr1, Chr = chr1, Freq = freq1,
                      n.snps.raw = n_snp_chr1, 
                      n.snps.refine = n_snp_chr1, 
                      cor.raw = 0, 
                      snp.start.raw = snp_start, snp.posStart.raw = snp.pos.start,
                      snp.end.raw = snp_end, snp.posEnd.raw = snp.pos.end,
                      cor.refine = 0, 
                      snp.start.refine = snp_start, snp.posStart.refine = snp.pos.start,
                      snp.end.refine = snp_end, snp.posEnd.refine = snp.pos.end,
                      flag.refine = "raw",
                      type.overlap.based.on.raw = "same",
                      stringsAsFactors = FALSE)
    res = rbind(res, res1)
    next
  }
  
  ## round1
  n_extend_round1 = 2*n_snp_chr1
  n_total = 5*n_snp_chr1
  if (n_total < 50) {
    n_extend_round1 = 25
  }
  
  idx_start_round1 = ifelse((idx_start - n_extend_round1) < 0, 1, (idx_start - n_extend_round1))
  flag_start_round1 = ifelse((idx_start - n_extend_round1) < 0, "out", "in")
  idx_end_round1 = ifelse((idx_end + n_extend_round1) > n_snps_cnvr1, n_snps_cnvr1, (idx_end + n_extend_round1))
  flag_end_round1 = ifelse((idx_end + n_extend_round1) > n_snps_cnvr1, "out", "in")
  
  
  ## set _chr1 for each cnvr1
  n_snp_round1 = idx_end_round1 - idx_start_round1 + 1
  mat_round1 = mat_LRR_cnvr1[, idx_start_round1:idx_end_round1]
  
  infor_max_round1 = refine_step1(Yt = mat_round1, min_len = 5)
  max.value.chr1 = infor_max_round1$max.value
  max.l.chr1 = infor_max_round1$max.l
  max.r.chr1 = infor_max_round1$max.r
  
  snp.start.round1 = colnames(mat_round1)[max.l.chr1]
  snp.end.round1   = colnames(mat_round1)[max.r.chr1]
  
  snp.round1.left  = colnames(mat_round1)[1]
  snp.round1.right = colnames(mat_round1)[n_snp_round1]
  
  ## round 2
  flag.round2 = NULL
  
  mat_round2 = NULL
  snp.start.round2 = NULL
  snp.end.round2 = NULL
  snp.round2.left = NULL
  snp.round2.right = NULL
  flag_start_round2 = NULL
  flag_end_round2 = NULL
  
  
  if ( (max.l.chr1 != 1) & (max.r.chr1 != n_snp_round1) ) {
    # donot do round2
    flag.round2 = "only_round1"
    
    mat_round2 = mat_round1
    n_snp_round2 = n_snp_round1
    snp.start.round2 = snp.start.round1
    snp.end.round2 = snp.end.round1
    
    snp.round2.left = snp.round1.left
    snp.round2.right = snp.round1.right
    
    flag_start_round2 = flag_start_round1
    flag_end_round2 = flag_end_round1
    
  } else if ( (max.l.chr1 == 1) & (max.r.chr1 != n_snp_round1) ) {
    flag.round2 = "leftside_round2"
    
    n_extend_round2 = 2*(n_extend_round1 + n_snp_chr1)
    idx_start_round2 = ifelse((idx_start - n_extend_round2) < 0, 1, (idx_start - n_extend_round2))
    flag_start_round2 = ifelse((idx_start - n_extend_round2) < 0, "out", "in")
    idx_end_round2 = idx_end_round1
    flag_end_round2 = flag_end_round1
    
    n_snp_round2 = idx_end_round2 - idx_start_round2 + 1
    mat_round2 = mat_LRR_cnvr1[, idx_start_round2:idx_end_round2]
    
    infor_max_round2 = refine_step1(Yt = mat_round2, min_len = 5)
    max.value.chr1 = infor_max_round2$max.value
    max.l.chr1 = infor_max_round2$max.l
    max.r.chr1 = infor_max_round2$max.r
    
    snp.start.round2 = colnames(mat_round2)[max.l.chr1]
    snp.end.round2 = colnames(mat_round2)[max.r.chr1]
    
    snp.round2.left = colnames(mat_round2)[1]
    snp.round2.right = colnames(mat_round2)[n_snp_round2]
    
  } else if ( (max.l.chr1 != 1) & (max.r.chr1 == n_snp_round1) ) {
    flag.round2 = "rightside_round2"
    
    n_extend_round2 = 2*(n_extend_round1 + n_snp_chr1)
    idx_start_round2 = idx_start_round1
    flag_start_round2 = flag_start_round1
    idx_end_round2 = ifelse((idx_end + n_extend_round2) > n_snp_cnvr1, n_snp_cnvr1, (idx_end + n_extend_round2))
    flag_end_round2 = ifelse((idx_end + n_extend_round2) > n_snp_cnvr1, "out", "in")
    
    n_snp_round2 = idx_end_round2 - idx_start_round2 + 1
    mat_round2 = mat_LRR_cnvr1[, idx_start_round2:idx_end_round2]
    
    infor_max_round2 = refine_step1(Yt = mat_round2, min_len = 5)
    max.value.chr1 = infor_max_round2$max.value
    max.l.chr1 = infor_max_round2$max.l
    max.r.chr1 = infor_max_round2$max.r
    
    snp.start.round2 = colnames(mat_round2)[max.l.chr1]
    snp.end.round2 = colnames(mat_round2)[max.r.chr1]
    
    snp.round2.left = colnames(mat_round2)[1]
    snp.round2.right = colnames(mat_round2)[n_snp_round2]
    
  } else if ( (max.l.chr1 == 1) & (max.r.chr1 == n_snp_round1) ) {
    flag.round2 = "bothside_round2"
    
    n_extend_round2 = 2*(n_extend_round1 + n_snp_chr1)
    idx_start_round2 = ifelse((idx_start - n_extend_round2) < 0, 1, (idx_start - n_extend_round2))
    flag_start_round2 = ifelse((idx_start - n_extend_round2) < 0, "out", "in")
    idx_end_round2 = ifelse((idx_end + n_extend_round2) > n_snp_cnvr1, n_snp_cnvr1, (idx_end + n_extend_round2))
    flag_end_round2 = ifelse((idx_end + n_extend_round2) > n_snp_cnvr1, "out", "in")
    
    n_snp_round2 = idx_end_round2 - idx_start_round2 + 1
    mat_round2 = mat_LRR_cnvr1[, idx_start_round2:idx_end_round2]
    
    infor_max_round2 = refine_step1(Yt = mat_round2, min_len = 5)
    max.value.chr1 = infor_max_round2$max.value
    max.l.chr1 = infor_max_round2$max.l
    max.r.chr1 = infor_max_round2$max.r
    
    snp.start.round2 = colnames(mat_round2)[max.l.chr1]
    snp.end.round2 = colnames(mat_round2)[max.r.chr1]
    
    snp.round2.left = colnames(mat_round2)[1]
    snp.round2.right = colnames(mat_round2)[n_snp_round2]
    
  }
  
  # results -----------------------------------------------------------------
  ## create annotate for heatmap
  idx1 = which(colnames(mat_round2) == snp_start)
  idx2 = which(colnames(mat_round2) == snp_end)
  # cnvr boundary refinement
  snps_flag_raw = c(rep("out_raw", idx1 - 1), 
                    rep("in_raw", (idx2-idx1+1)), 
                    rep("out_raw", (n_snp_round2 - idx2)))
  stopifnot(max.l.chr1 < max.r.chr1) ## check index
  ## round1 
  idx.l.round1 = which(colnames(mat_round2) == snp.start.round1)
  idx.r.round2 = which(colnames(mat_round2) == snp.end.round1)
  snps_flag_refine_round1 = c(rep("out_refine_round1", idx.l.round1 - 1), 
                              rep("in_refine_round1", (idx.r.round2 - idx.l.round1 + 1)), 
                              rep("out_refine_round1", (n_snp_round2 - idx.r.round2)))
  
  ## round2
  snps_flag_refine_round2 = c(rep("out_refine_round2", max.l.chr1 - 1),
                              rep("in_refine_round2", (max.r.chr1 - max.l.chr1 + 1)),
                              rep("out_refine_round2", (n_snp_round2 - max.r.chr1)))
  
  ## final result
  snps_flag_final = NULL
  idx.start.final = max.l.chr1
  idx.end.final = max.r.chr1
  if (max.l.chr1 == 1 | max.r.chr1 == n_snp_round2) { 
    ## can be changed here if oneside in the round2 to the left or right 
    ## use raw CNVR boundary
    snps_flag_final = c(rep("out_final", idx1 - 1), 
                        rep("in_final", (idx2-idx1+1)), 
                        rep("out_final", (n_snp_round2 - idx2)))
    idx.start.final = idx1
    idx.end.final = idx2
  } else {
    snps_flag_final = c(rep("out_final", max.l.chr1 - 1),
                        rep("in_final", (max.r.chr1 - max.l.chr1 + 1)),
                        rep("out_final", (n_snp_round2 - max.r.chr1)))
  }
  
  ## annotate in round1 and round2
  snps_final = colnames(mat_round2)
  n.left.round2 = which(snps_final == snp.round1.left) - which(snps_final == snp.round2.left)
  n.right.round2 = which(snps_final == snp.round2.right) - which(snps_final == snp.round1.right)
  
  n.left.round1 = which(snps_final == snp_start) - which(snps_final == snp.round1.left)
  n.right.round1 = which(snps_final == snp.round1.right) - which(snps_final == snp_end)
  
  n.raw = n_snp_round2 - n.left.round1 - n.left.round2 - n.right.round1 - n.right.round2
  snp_flag_round = c(rep("round2_add", n.left.round2),
                     rep("round1_add", n.left.round1),
                     rep("raw", n.raw),
                     rep("round1_add", n.right.round1),
                     rep("round2_add", n.right.round2))
  
  ## add information
  mcorr_round2 = cor(mat_round2)
  n2 = length(mcorr_round2[lower.tri(mcorr_round2)])
  S = mcorr_round2[lower.tri(mcorr_round2)]
  
  mcorr_raw = cor( mat_round2[, idx1:idx2])
  x = mcorr_raw[lower.tri(mcorr_raw)]
  n1 = length(x)
  max.value.raw= (mean(x) - (sum(S)-sum(x))/(n2-n1))/sqrt(1/n1+1/(n2-n1))
  
  if (flag.round2 == "only_round1") {
    annotate_col1 = data.frame(
      group_raw = snps_flag_raw,
      group_refine_round1 = snps_flag_refine_round1,
      group_round_add = snp_flag_round,
      group_final = snps_flag_final
    )
  } else {
    annotate_col1 = data.frame(
      group_raw = snps_flag_raw,
      group_refine_round1 = snps_flag_refine_round1,
      group_refine_round2 = snps_flag_refine_round2,
      group_round_add = snp_flag_round,
      group_final = snps_flag_final
    )
  }
  
  rownames(annotate_col1) = colnames(mat_round2)
  
  snp_cnvr1_round2 = colnames(mat_round2) # for following res1
  
  ## add snp_position information
  # for raw cnvr1 boundary
  snp.start.raw.cnvr1 = snp_cnvr1_round2[idx1]
  snp.end.raw.cnvr1 = snp_cnvr1_round2[idx2]
  snp.posStart.raw.cnvr1 = snp_chr1$Position[which(snp_chr1$Name == snp.start.raw.cnvr1)]
  snp.posEnd.raw.cnvr1 = snp_chr1$Position[which(snp_chr1$Name == snp.end.raw.cnvr1)]
  snp.num.raw.cnvr1 = idx2 - idx1 + 1
  # for refine cnvr1 boundary
  snp.start.refine.cnvr1 = snp_cnvr1_round2[idx.start.final]
  snp.end.refine.cnvr1 = snp_cnvr1_round2[idx.end.final]
  snp.posStart.refine.cnvr1 = snp_chr1$Position[which(snp_chr1$Name == snp.start.refine.cnvr1)]
  snp.posEnd.refine.cnvr1 = snp_chr1$Position[which(snp_chr1$Name == snp.end.refine.cnvr1)]
  snp.num.refine.cnvr1 = idx.end.final - idx.start.final + 1
  
  # check overlaptype.refine.raw 
  type.overlap.cnvr1 = NULL
  if (idx1 == idx.start.final & idx2 == idx.end.final) {
    type.overlap.cnvr1 = "same"
  } else if (idx.end.final <= idx1) {
    type.overlap.cnvr1 = "out_left"
  } else if (idx.start.final >= idx2) {
    type.overlap.cnvr1 = "out_right"
  } else if ((idx.start.final <= idx1) & (idx.end.final >= idx2)) {
    type.overlap.cnvr1 = "include_out"
  } else if ((idx1 <= idx.start.final) & (idx2 >= idx.end.final)) {
    type.overlap.cnvr1 = "include_in"
  } else if (idx.start.final<idx1 & idx.end.final>idx1 & idx.end.final<= idx2) {
    type.overlap.cnvr1 = "overlvep_left"
  } else if (idx.end.final>idx2 & idx.start.final>=idx1 & idx.start.final<idx2) {
    type.overlap.cnvr1 = "overlap_right"
  }
  
  # pheatmap
  if ( flag_plot ) {
    filename_png = paste0(cnvr1, "_boundary_refinement.png")
    png(filename = file.path(path_png, filename_png), width = 12, height = 12, res = 512, units = "in")
    par(mar = c(4, 4, 4, 4))
    pheatmap(mcorr_round2,
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             annotation_col = annotate_col1,
             show_rownames = FALSE,
             show_colnames = FALSE,
             main = paste(cnvr1, "type_refine:", flag.round2, 
                          "flag_boundary:", paste0("(left:", flag_start_round2,")", "(right:", flag_end_round2, ")"),
                          "\n", "type.overlap:", type.overlap.cnvr1,
                          "freq:", freq1,
                          "cor.raw:", round(max.value.raw, 2), 
                          "cor.refine:", round(max.value.chr1, 2)))
    dev.off()
  }

  res1 = data.frame(CNVR_ID = cnvr1, 
                    Chr = chr1, 
                    Freq = freq1, 
                    n.snps.raw = snp.num.raw.cnvr1, 
                    n.snps.refine = snp.num.refine.cnvr1, 
                    cor.raw = max.value.raw, 
                    snp.start.raw = snp.start.raw.cnvr1, 
                    snp.posStart.raw = snp.posStart.raw.cnvr1,
                    snp.end.raw = snp.end.raw.cnvr1, 
                    snp.posEnd.raw = snp.posEnd.raw.cnvr1,
                    cor.refine = max.value.chr1, 
                    snp.start.refine = snp.start.refine.cnvr1, 
                    snp.posStart.refine = snp.posStart.refine.cnvr1,
                    snp.end.refine = snp.end.refine.cnvr1, 
                    snp.posEnd.refine = snp.posEnd.refine.cnvr1,
                    flag.refine = "refine",
                    type.overlap.based.on.raw = type.overlap.cnvr1,
                    stringsAsFactors = FALSE)
  res = rbind(res, res1)
  
} 

filename = paste0("CNVR_refine_chr_", chr1, "_detail.rds")
saveRDS(res, file = file.path(path_res, filename))




