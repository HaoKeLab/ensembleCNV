#!/hpc/packages/minerva-common/R/3.3.1/lib64/R/bin/Rscript --vanilla

library(data.table)
library(tibble)

file_sampleID_transform = "/sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/IPQ/IPQ_batch2_split/dat/sampleID_transform.rds"
dat_sampleID_transform = readRDS(file = file_sampleID_transform)

read_matrix = function(path_input, filename_matrix, dat_sampleID_transform, b1) {
  
  dat_LRR = fread(input = file.path( path_input, filename_matrix) )
  
  dat_LRR <- as.data.frame(dat_LRR, stringsAsFactors = FALSE)
  dat_LRR <- column_to_rownames(dat_LRR, var = "V1")
  
  sampleID_raw = rownames( dat_LRR )
  dat1 = subset(dat_sampleID_transform, batch == b1)
  
  sampleID_new = dat1$Sample_ID_new[match(sampleID_raw, dat1$Sample_ID_file)]
  
  stopifnot( length(na.omit(sampleID_new)) == length(sampleID_raw) )
  
  rownames(dat_LRR) = sampleID_new
  
  dat_LRR
}


path_input = "/sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/SNP_recluster/dat"

dat_batch1 = read_matrix(path_input = path_input, filename_matrix = "matrix.LRR.batch1.txt", dat_sampleID_transform = dat_sampleID_transform, b1 = 1)

dat_batch2_1 = read_matrix(path_input = path_input, filename_matrix = "matrix.LRR.batch2.1.txt", dat_sampleID_transform = dat_sampleID_transform, b1 = "2_1")
dat_batch2_2 = read_matrix(path_input = path_input, filename_matrix = "matrix.LRR.batch2.2.txt", dat_sampleID_transform = dat_sampleID_transform, b1 = "2_2")
dat_batch2_3 = read_matrix(path_input = path_input, filename_matrix = "matrix.LRR.batch2.3.txt", dat_sampleID_transform = dat_sampleID_transform, b1 = "2_3")

dat_batch3 = read_matrix(path_input = path_input, filename_matrix = "matrix.LRR.batch3.txt", dat_sampleID_transform = dat_sampleID_transform, b1 = 3)


dat_all = rbind(dat_batch1, dat_batch2_1,
                dat_batch2_2, dat_batch2_3,
                dat_batch3)

sampleID = rownames(dat_all)
stopifnot( length(unique(sampleID)) == nrow(dat_all) )  ##

mat = as.matrix(dat_all)
rownames(mat) = sampleID
colnames(mat) = NULL

dim(mat)
mat[1:5, 1:5]
## output path 
path_output = "/sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/SNP_recluster/"
saveRDS(mat, file = file.path(path_output, "matrix.LRR.all.rds"))

col_mean = colMeans(mat, na.rm = TRUE)
for (i in 1:nrow(mat)) {
  cat(i, "\n")
  v1 = as.vector(mat[i, ])
  idx1 = which(is.na(v1))
  if (length(idx1) >= 1) {
    mat[i, idx1] = col_mean[idx1]
  }
}

dat = as.data.frame(mat)
rownames(dat) = sampleID
dat[1:5, 1:5]

PCA = prcomp(dat)
PC  = predict(PCA)
PC  = as.data.frame(PC, stringsAsFactors = FALSE)
PC  = PC[, c("PC1", "PC2", "PC3")]

path_output = "/sc/orga/projects/haok01a/chengh04/paper/ensembleCNV/Fig_batch_effect_FA/data"
saveRDS(PC, file = file.path(path_output, "PC_snps_after_QC.png"))




























# path_LRR <- "/sc/orga/projects/haok01a/chengh04/MEGA/MegaEX_Inga/analysis_Part1/batch_effect/res"
# filename_LRR <- "matrix.LRR.snps.randomly.select.txt"
# dat_LRR <- fread(input = file.path(path_LRR, filename_LRR))
# 
# dat_LRR <- as.data.frame(dat_LRR, stringsAsFactors = FALSE)
# dat_LRR <- column_to_rownames(dat_LRR, var = "V1")
# 
# mat_LRR <- as.matrix(dat_LRR)
# colnames(mat_LRR) <- NULL
# 
# ## fillin NA value
# median_LRR <- unlist(lapply(1:ncol(mat_LRR), FUN = function(k) {
#   median(mat_LRR[, k], na.rm = TRUE)
# }))
# 
# for (i in 1:nrow(mat_LRR)) {
#   v1 <- mat_LRR[i, ]
#   idx.na <- which(is.na(v1))
#   if (length(idx.na) >= 1) {
#     mat_LRR[i, idx.na] <- median_LRR[idx.na]
#   }
# }
# 
# ## PCA
# PCA <- prcomp(x = mat_LRR)
# PC  <- predict(PCA)
# 
# saveRDS(PCA, file = file.path(path_LRR, "PCA.rds"))
# saveRDS(PC, file = file.path(path_LRR, "PC.rds"))
# 
# 
# ###################################################################
# anno.samples <- read.table(file = "/sc/orga/projects/haok01a/chengh04/MEGA/MegaEX_Inga/FinalReport/Plates_Part1/MegaEX_Inga_Part1_Samples_Table.txt",
#                            header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE, na.strings = "")
# samples.PC   <- rownames(PC)
# samples.anno <- anno.samples$`Sample ID`
# length( samples.PC )
# length( samples.anno )
# stopifnot( length(samples.PC) == length(samples.anno) )
# sum( samples.PC %in% samples.anno)
# 
# ## -------------------------------------------------------
# # idx1 <- which(samples.anno %in% samples.PC)
# # dat1 <- anno.samples[idx1, ]
# # dat1$Uniq_Sample_ID <- dat1$`Sample ID`
# # nrow(dat1)
# # dat2 <- anno.samples[-idx1, ]
# # dat2$Uniq_Sample_ID <- paste(dat2$`Sample ID`, 
# #                              dat2$`Array Info.Sentrix ID`,
# #                              dat2$`Array Info.Sentrix Position`, sep = "_")
# # nrow(dat2)
# # sum(dat2$Uniq_Sample_ID %in% samples.PC)
# # 
# # anno.samples <- rbind(dat1, dat2)
# ###############################################
# anno.samples <- anno.samples[, c("Sample ID", "Sample Plate")]
# names(anno.samples) <- c("Sample_ID", "Sample_Plate")
# anno.samples$Sample_Plate[which(is.na(anno.samples$Sample_Plate))] <- "Plate34&35"
# 
# PC.annotate <- merge(PC, anno.samples,  ## 1859 
#                      by.x = "row.names", by.y = "Sample_ID")
# saveRDS(PC.annotate, file = file.path(path_LRR, "PC.annotate.rds"))
# 
# library(gridExtra)
# library(ggplot2)
# 
# ## PCA
# # png(filename = file.path(path_LRR, "PCA.png"), width = 12, height = 12, units = "in", res = 512)
# # plot(PCA)
# # dev.off()
# 
# PCA.summary <- summary(PCA)
# PCA.important <- PCA.summary$importance
# PCA.important.pc <- PCA.important[3, 1:10]
# dat.cum.prop <- data.frame(PC = 1:length(PCA.important.pc), 
#                            prop = PCA.important.pc)
# rownames(dat.cum.prop) <- NULL
# 
# ## Cumulative Proportion 
# png(filename = file.path(path_LRR, "PCA_importance.png"), width = 12, height = 12, units = "in", res = 512)
# g1 <- ggplot(data = dat.cum.prop, aes(PC, prop, group = 1)) + 
#   geom_point(size = 4) + 
#   geom_line() + 
#   theme_bw(base_size = 9) +
#   theme(plot.title = element_text(size = 15)) +
#   ggtitle("Cumulative Proportion for first 10 PCs")
# print(g1)
# dev.off()
# 
# png(filename = file.path(path_LRR, "PCA_PC12_annotate.png"),
#     width = 12, height = 12, units = "in", res = 512)
# p1 <- ggplot(data = PC.annotate, aes(PC1, PC2, col = Sample_Plate)) + 
#   geom_point(size = 3) + 
#   theme_bw() + 
#   theme(legend.title = element_text(size = 13),
#         legend.text = element_text(size = 13),
#         plot.title = element_text(size = 15),
#         legend.position = "top") + 
#   ggtitle("PC2 ~ PC1")
# print(p1)
# dev.off()
# 
# png(filename = file.path(path_LRR, "PCA_PC13_annotate.png"),
#     width = 12, height = 12, units = "in", res = 512)
# p2 <- ggplot(data = PC.annotate, aes(PC1, PC3, col = Sample_Plate)) + 
#   geom_point(size = 3) + 
#   theme_bw() + 
#   theme(legend.title = element_text(size = 13),
#         legend.text = element_text(size = 13),
#         plot.title = element_text(size = 15),
#         legend.position = "top") + 
#   ggtitle("PC3 ~ PC1")
# print(p2)
# dev.off()
# 
# png(filename = file.path(path_LRR, "PCA_PC23_annotate.png"),
#     width = 12, height = 12, units = "in", res = 512)
# p3 <- ggplot(data = PC.annotate, aes(PC2, PC3, col = Sample_Plate)) + 
#   geom_point(size = 3) + 
#   theme_bw() + 
#   theme(legend.title = element_text(size = 13),
#         legend.text = element_text(size = 13),
#         plot.title = element_text(size = 15),
#         legend.position = "top") + 
#   ggtitle("PC3 ~ PC2")
# print(p3)
# dev.off()
# 
# 
# # for clear only annotate Plates 75 & 76 ----------------------------------
# 
# table(PC.annotate$Sample_Plate) # `Plate 75` & `Plate 76`
# 
# PC.annotate.clear <- PC.annotate
# PC.annotate.clear$Plate_annotate <- PC.annotate.clear$Sample_Plate
# PC.annotate.clear$Plate_annotate[which(!PC.annotate.clear$Plate_annotate %in%  c("Plate 75", "Plate 76"))] <- "Plate_other"
# table(PC.annotate.clear$Plate_annotate) # `Plate 75` & `Plate 76`
# 
# 
# png(filename = file.path(path_LRR, "PCA_PC12_annotate_clear.png"),
#     width = 12, height = 12, units = "in", res = 512)
# p1 <- ggplot(data = PC.annotate.clear, aes(PC1, PC2, col = Plate_annotate)) + 
#   geom_point(size = 3) + 
#   theme_bw() + 
#   theme(legend.title = element_text(size = 13),
#         legend.text = element_text(size = 13),
#         plot.title = element_text(size = 15),
#         legend.position = "top") + 
#   ggtitle("PC2 ~ PC1")
# print(p1)
# dev.off()
# 
# png(filename = file.path(path_LRR, "PCA_PC13_annotate_clear.png"),
#     width = 12, height = 12, units = "in", res = 512)
# p2 <- ggplot(data = PC.annotate.clear, aes(PC1, PC3, col = Plate_annotate)) + 
#   geom_point(size = 3) + 
#   theme_bw() + 
#   theme(legend.title = element_text(size = 13),
#         legend.text = element_text(size = 13),
#         plot.title = element_text(size = 15),
#         legend.position = "top") + 
#   ggtitle("PC3 ~ PC1")
# print(p2)
# dev.off()
# 
# png(filename = file.path(path_LRR, "PCA_PC23_annotate_clear.png"),
#     width = 12, height = 12, units = "in", res = 512)
# p3 <- ggplot(data = PC.annotate.clear, aes(PC2, PC3, col = Plate_annotate)) + 
#   geom_point(size = 3) + 
#   theme_bw() + 
#   theme(legend.title = element_text(size = 13),
#         legend.text = element_text(size = 13),
#         plot.title = element_text(size = 15),
#         legend.position = "top") + 
#   ggtitle("PC3 ~ PC2")
# print(p3)
# dev.off()
# 
# 
