## 2016-10-07 this is the second version of PCA method
# add all the supplement in this code

rm(list = ls())

library(ggplot2)
library(gdata)
library(gridExtra)
## load data
samples.file <- "E:/Desktop/Food_Allergy/cnv-pca-qc/data/QC.samples.new.rds"
plink.pca.file <- "E:/Desktop/Food_Allergy/cnv-pca-qc/data/plink.pca.tab"
batch.cnv.file = "E:/Desktop/Food_Allergy/cnv-pca-qc/data/batch for cnv.xls"

batch.cnv <- read.xls(batch.cnv.file)  ## batch samples data
batch.cnv$iid <- as.character(batch.cnv$iid)
batch.cnv$duplication <- as.character(batch.cnv$duplication)
batch.cnv$SentrixPosition_A <- as.character(batch.cnv$SentrixPosition_A)
idxs1 <- which(!batch.cnv$duplication %in% c("database", "duplicate"))
batch.cnv$Sample_ID[idxs1] = batch.cnv$iid[idxs1]
idxs2 <- which(batch.cnv$duplication %in% c("database", "duplicate"))
batch.cnv$Sample_ID[idxs2] = paste(batch.cnv$iid[idxs2], 
                                   batch.cnv$SentrixBarcode_A[idxs2],
                                   batch.cnv$SentrixPosition_A[idxs2], sep = "_")

dt <- readRDS(file = samples.file) ## samples data
dt.plink.pca <- read.table(file = plink.pca.file,  ## antonio' pca data
                     header = TRUE, sep = "\t",
                     colClasses = c("NULL", "character", rep("numeric", 10)))

## batch data 
file_batch <- "C:/HaoxiangCheng/Food.Allergy/FA/code_batch/analysis_batch/dat/samples.anno.batch.gs.txt"
dat1  <- read.table(file = file_batch,
                   sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
dat1$Sample_ID <- gsub(":", "_", dat1$Sample.Uniq.ID)  ## for 3 batch compare
dat1$batch <- ifelse(dat1$batch.recluster == 1, "batch1",
                     ifelse(dat1$batch.recluster == 2, "batch2", "batch3"))
dat_batch  <- dat1[, c("Sample_ID", "batch")]
table(dat_batch$batch)

samples_batch1 <- subset(dat_batch, batch == "batch1")$Sample_ID
samples_batch2 <- subset(dat_batch, batch == "batch2")$Sample_ID
samples_batch3 <- subset(dat_batch, batch == "batch3")$Sample_ID


## ----------------------------------------------------------------------
# dt.normal
features.normal <- function(dt) {
  dt <- lapply(dt, FUN = function(x) (x - median(x))/mad(x))
  dtt <- do.call(cbind, dt)
  dtt
}

# PAC
run.pca <- function(dt, features, tt) {
  
  dt1 <- dt[, features]
  dt1 <- features.normal(dt = dt1)
  
  pca <- princomp(dt1) 
  
  pdf(file = file.path('C:/HaoxiangCheng/Food.Allergy/FA/code_batch/analysis_batch/png/', paste0(tt, "-time-pca-part1.pdf")))
  
  plot(pca, main = paste(tt, "PAC"))
  biplot(pca, main = paste(tt, "PAC"))
  
  pca.result <- data.frame(pca$scores, stringsAsFactors = FALSE)
  
  p1 = ggplot(pca.result, aes(Comp.1, Comp.2)) +
    geom_point() 

  p2 = ggplot(pca.result, aes(Comp.1, Comp.3)) +
    geom_point()

  p3 = ggplot(pca.result, aes(Comp.2, Comp.3)) +
    geom_point() 

  p4 = ggplot(pca.result, aes(Comp.2, Comp.3)) +
    geom_blank() +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          line = element_blank(),
          panel.background = element_blank())

  gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)
  
  dev.off()
  dt.out <- cbind(dt, pca.result)
  dt.out
} 


summary.pca.results.first <- function(dt, comp1.cutoff, comp2.cutoff = NULL, tt, dt.plink.pca, batch.cnv) {
  
  pdf(file = file.path('C:/HaoxiangCheng/Food.Allergy/FA/code_batch/analysis_batch/png/', paste0(tt, "-time-pca-part2.pdf")))
  
  plot(dt$Comp.1, dt$Comp.2, main = paste(tt, "time PCA"))
  abline(v = comp1.cutoff, h = comp2.cutoff)

  # dat_plate <- batch.cnv[, c("Sample_Plate", "Sample_ID")]
  # dat_plate$Plate <- ifelse(dat_plate$Sample_Plate == "XW_FA_P30", "XW_FA_P30",
  #                           ifelse(dat_plate$Sample_Plate == "XW_FA_P31", "XW_FA_P31",
  #                                  ifelse(dat_plate$Sample_Plate == "XW_FA_P32", "XW_FA_P32",
  #                                         ifelse(dat_plate$Sample_Plate == "XW_FA_P33", "XW_FA_P33", "other"))))
  # dt1 <- merge(dt, dat_plate)
  # 
  # g1 <- ggplot(dt1, aes(Comp.1, Comp.2, col = Plate)) + 
  #   geom_point(shape = 2) + 
  #   geom_hline(yintercept = comp2.cutoff, lty = 2) + 
  #   geom_vline(xintercept = comp1.cutoff, lty = 2) 
  # 
  # print(g1)
  
  idx1 <- which(dt$Comp.1 < comp1.cutoff)
  if (is.null(comp2.cutoff)) {
    idx2 <- NULL
  }  else {
    idx2 <- which(dt$Comp.2 < comp2.cutoff)
  }
  
  
  idxs <- union(idx1, idx2)
  
  dt.filter <- dt[idxs, ]
  dt.keep <- dt[-idxs, ]
  
  sink(file = file.path('C:/HaoxiangCheng/Food.Allergy/FA/code_batch/analysis_batch/png/', paste0(tt, "-time-pca.txt")))
  
  ## summary of PennCNV.BAF_drift
  cat("samples-filter-BAF-drift-summary:\n")
  print(summary(dt.filter$PennCNV.BAF_drift))
  cat("samples-keep-BAF-drift-summary:\n")
  print(summary(dt.keep$PennCNV.BAF_drift))
  
  cat("samples-filter-PennCNV-NumCNV:\n")
  print(summary(dt.filter$PennCNV.NumCNV))
  cat("samples-keep-PennCNV-NumCNV:\n")
  print(summary(dt.keep$PennCNV.NumCNV))
  cat("samples-filter-QuantiSNP-NumCNV:\n")
  print(summary(dt.filter$QuantiSNP.NumCNV))
  cat("samples-keep-QuantiSNP-NumCNV:\n")
  print(summary(dt.keep$QuantiSNP.NumCNV))
  cat("samples-filter-iPattern-NumCNV:\n")
  print(summary(dt.filter$iPattern.NumCNV))
  cat("samples-keep-iPattern-NumCNV:\n")
  print(summary(dt.keep$iPattern.NumCNV))
  
  samples.filter <- dt.filter$IID
  samples.keep <- dt.keep$IID
  
  cat("number of filter samples: \n")
  print(length(samples.filter))
  cat("number of keep samples: \n")
  print(length(samples.keep))
  
  ## family information
  cat("families of samples-filter:\n")
  print(table(substr(samples.filter, 1, 6)))
  
  ## plot above samples distribution of antonio's pc results
  dt.plink.pca$cluster = "1"
  idxs = which(dt.plink.pca$IID %in% samples.filter)
  dt.plink.pca$cluster[idxs] = "2"
  dt.pca1 = subset(dt.plink.pca, cluster == 1)
  dt.pca2 = subset(dt.plink.pca, cluster == 2)
  
  p5 = ggplot(data = dt.pca1, aes(PC1, PC2, colour = 'other')) + 
    geom_point() +
    geom_point(data = dt.pca2, aes(PC1, PC2, colour = "filter")) + 
    ggtitle(paste(tt, "time PAC data distribution"))
  gridExtra::grid.arrange(p5)
  
  ## 
  idxs = which(batch.cnv$Sample_ID %in% samples.filter)
  batch.cnv.sub = batch.cnv[idxs, ]
  
  cat("samples.filter Plate:\n")
  print(table(batch.cnv.sub$plate))
  cat("All Samples Plate:\n")
  print(table(batch.cnv$plate))
  
  dev.off()
  sink()
  return(list(samples.filter = samples.filter,
              samples.keep = samples.keep))
  
}

summary.pca.results.second <- function(dt, comp1.cutoff, comp2.cutoff = NULL, tt, dt.plink.pca, batch.cnv,
                                       samples_batch1, samples_batch2) {
  
  pdf(file = file.path('C:/HaoxiangCheng/Food.Allergy/FA/code_batch/analysis_batch/png/', paste0(tt, "-time-pca-part2.pdf")))
  
  # plot(dt$Comp.1, dt$Comp.2, main = paste(tt, "time PCA"))
  # abline(v = comp1.cutoff, h = comp2.cutoff)
  # 
  dt1 <- dt[, c("IID","Comp.1", "Comp.2")]
  names(dt1)[1] <- "Sample_ID"
  dat_plate <- batch.cnv[, c("Sample_Plate", "Sample_ID")]
  dat_plate$Plate <- ifelse(dat_plate$Sample_Plate == "XW_FA_P30", "XW_FA_P30",
                            ifelse(dat_plate$Sample_Plate == "XW_FA_P31", "XW_FA_P31",
                                   ifelse(dat_plate$Sample_Plate == "XW_FA_P32", "XW_FA_P32",
                                          ifelse(dat_plate$Sample_Plate == "XW_FA_P33", "XW_FA_P33", "other"))))
  
  dat_annotate = readRDS(file = "E:/Desktop/Food_Allergy/cnv-pca-qc/FA_2790_annotate.rds")
  
  dt1 <- merge(dt1, dat_annotate)
  nrow(dt1)

  col = c("KEEP" = "gray", "XW-FA-P04" = "red", 
          "XW_FA_P30" = "green3", "XW_FA_P31" = "blue",
          "XW_FA_P32" = "cyan", "XW_FA_P33" = "magenta",
          "OUTLIER" = "black")
  
  dt1$flag_Sample = dt1$flag
  
  png(filename = file.path('C:/HaoxiangCheng/Food.Allergy/FA/code_batch/analysis_batch/png/second_time_PCA.png'),
      width = 12, height = 12, res = 512, units = "in")
  g1 <- ggplot(dt1, aes(Comp.1, Comp.2, col = flag_Sample)) + 
    geom_point(shape = 1, size = 3) + 
    # geom_vline(xintercept = comp1.cutoff, lty = 2) +
    theme_bw(base_size = 9) +
    theme(legend.position = "top",
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 15)) +
    xlab("PC1") + 
    ylab("PC2") + 
    scale_color_manual(values = col)
  
  print(g1)
  
  dev.off()
  
  dt1$Plate[which(dt1$Sample_ID %in% samples_batch1)] <- "batch1"
  dt1$Plate[which(dt1$Sample_ID %in% samples_batch2)] <- "batch2"
  
  g2 <- ggplot(dt1, aes(Comp.1, Comp.2, col = Plate)) + 
    geom_point(shape = 2) + 
    geom_vline(xintercept = comp1.cutoff, lty = 2) +
    theme_bw(base_size = 9) +
    theme(legend.position = "top")
  
  print(g2)
  
  idx1 <- which(dt$Comp.1 < comp1.cutoff)
  if (is.null(comp2.cutoff)) {
    idx2 <- NULL
  }  else {
    idx2 <- which(dt$Comp.2 < comp2.cutoff)
  }
  
  
  idxs <- union(idx1, idx2)
  
  dt.filter <- dt[idxs, ]
  dt.keep <- dt[-idxs, ]
  
  sink(file = file.path('C:/HaoxiangCheng/Food.Allergy/FA/code_batch/analysis_batch/png/', paste0(tt, "-time-pca.txt")))

  ## summary of PennCNV.BAF_drift
  cat("samples-filter-BAF-drift-summary:\n")
  print(summary(dt.filter$PennCNV.BAF_drift))
  cat("samples-keep-BAF-drift-summary:\n")
  print(summary(dt.keep$PennCNV.BAF_drift))
  
  cat("samples-filter-PennCNV-NumCNV:\n")
  print(summary(dt.filter$PennCNV.NumCNV))
  cat("samples-keep-PennCNV-NumCNV:\n")
  print(summary(dt.keep$PennCNV.NumCNV))
  cat("samples-filter-QuantiSNP-NumCNV:\n")
  print(summary(dt.filter$QuantiSNP.NumCNV))
  cat("samples-keep-QuantiSNP-NumCNV:\n")
  print(summary(dt.keep$QuantiSNP.NumCNV))
  cat("samples-filter-iPattern-NumCNV:\n")
  print(summary(dt.filter$iPattern.NumCNV))
  cat("samples-keep-iPattern-NumCNV:\n")
  print(summary(dt.keep$iPattern.NumCNV))
  
  samples.filter <- dt.filter$IID
  samples.keep <- dt.keep$IID
  
  cat("number of filter samples: \n")
  print(length(samples.filter))
  cat("number of keep samples: \n")
  print(length(samples.keep))
  
  ## family information
  cat("families of samples-filter:\n")
  print(table(substr(samples.filter, 1, 6)))
  
  ## plot above samples distribution of antonio's pc results
  dt.plink.pca$cluster = "1"
  idxs = which(dt.plink.pca$IID %in% samples.filter)
  dt.plink.pca$cluster[idxs] = "2"
  dt.pca1 = subset(dt.plink.pca, cluster == 1)
  dt.pca2 = subset(dt.plink.pca, cluster == 2)
  
  p5 = ggplot(data = dt.pca1, aes(PC1, PC2, colour = 'other')) + 
    geom_point() +
    geom_point(data = dt.pca2, aes(PC1, PC2, colour = "filter")) + 
    ggtitle(paste(tt, "time PAC data distribution"))
  gridExtra::grid.arrange(p5)
  
  ## 
  idxs = which(batch.cnv$Sample_ID %in% samples.filter)
  batch.cnv.sub = batch.cnv[idxs, ]
  
  cat("samples.filter Plate:\n")
  print(table(batch.cnv.sub$plate))
  cat("All Samples Plate:\n")
  print(table(batch.cnv$plate))
  
  dev.off()
  sink()
  return(list(samples.filter = samples.filter,
              samples.keep = samples.keep))
  
}



## -------------------------------------------------------------------------------------

## ----------------------first time pca
# features name for PCA
features <- c("PennCNV.LRR_SD", "PennCNV.BAF_SD", "PennCNV.BAF_drift", "PennCNV.WF",
              "PennCNV.NumCNV", "QuantiSNP.LRR_SD", "QuantiSNP.BAF_SD",
              "QuantiSNP.NumCNV", "iPattern.NumCNV", "iPattern.LRR_SD")

## plot PCA results 
pca.out <- run.pca(dt = dt, features = features, tt = "first")
## summary pca results
pca.res <- summary.pca.results.first(dt = pca.out, comp1.cutoff = -15, comp2.cutoff = -15, tt ="first",
                                     dt.plink.pca = dt.plink.pca, batch.cnv = batch.cnv)

saveRDS(pca.res$samples.filter, file = 'C:/HaoxiangCheng/Food.Allergy/FA/code_batch/analysis_batch/png/first-samples-filter.rds')
saveRDS(pca.res$samples.keep, file = 'C:/HaoxiangCheng/Food.Allergy/FA/code_batch/analysis_batch/png/first-samples-keep.rds')


## ---------------------second time pca
samples.filter <- pca.res$samples.filter
dt1 <- dt[which(!dt$IID %in% samples.filter), ]

# feature selection
features1 <- c("PennCNV.LRR_SD", "PennCNV.BAF_SD", "PennCNV.WF",
              "PennCNV.NumCNV", "QuantiSNP.LRR_SD", "QuantiSNP.BAF_SD",
              "QuantiSNP.NumCNV", "iPattern.NumCNV", "iPattern.LRR_SD")

## plot PCA results 
pca.out <- run.pca(dt = dt1, features = features1, tt = "second")
## summary pca results
pca.res <- summary.pca.results.second(dt = pca.out, comp1.cutoff = -3, tt ="second",
                                      dt.plink.pca = dt.plink.pca, batch.cnv = batch.cnv,
                                      samples_batch1 = samples_batch1,
                                      samples_batch2 = samples_batch2)

saveRDS(pca.res$samples.filter, file = 'C:/HaoxiangCheng/Food.Allergy/FA/code_batch/analysis_batch/png/second-samples-filter.rds')
saveRDS(pca.res$samples.keep, file = 'C:/HaoxiangCheng/Food.Allergy/FA/code_batch/analysis_batch/png/second-samples-keep.rds')




## ----------------------third time pca
# samples.filter <- pca.res$samples.filter
# dt2 <- dt1[which(!dt1$IID %in% samples.filter), ]
# 
# features2 <- c("PennCNV.LRR_SD", "PennCNV.BAF_SD", "PennCNV.WF",
#                 "QuantiSNP.LRR_SD", "QuantiSNP.BAF_SD",
#                  "iPattern.LRR_SD")
# 
# ## plot PCA results 
# pca.out <- run.pca(dt = dt2, features = features2, tt = "third")
# pca.res <- summary.pca.results(dt = pca.out, comp1.cutoff = -7, tt ="third",
#                                dt.plink.pca = dt.plink.pca, batch.cnv = batch.cnv)

