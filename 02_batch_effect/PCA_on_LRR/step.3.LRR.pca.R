#!/usr/bin/env Rscript

## PCA
args <- commandArgs( trailingOnly = TRUE )
wk_dir <- args[1]            ## working directory where the LRR matrix is located for PCA
filename_matrix <- args[2]   ## the LRR matrix generated in step 2

suppressMessages({
  require(data.table)
  require(tibble)
  require(cowplot)
  require(ggplot2)
})


  dat_LRR <- fread(input = file.path( wk_dir, filename_matrix) )  
  dat_LRR <- as.data.frame(dat_LRR, stringsAsFactors = FALSE)
  dat_LRR <- column_to_rownames(dat_LRR, var = "V1")
  
  sampleID <- rownames( dat_LRR )
  
  ## deal with NA values in matrix
  mat <- as.matrix(dat_LRR)
  rownames(mat) <- sampleID
  colnames(mat) <- NULL
  
  col_mean <- colMeans(mat, na.rm = TRUE)
  for (i in 1:nrow(mat)) {
    v1 <- as.vector(mat[i, ])
    idx1 <- which(is.na(v1))
    if (length(idx1) >= 1) {
      mat[i, idx1] <- col_mean[idx1]
    }
  }
  
  ## check which SNPs with all values being NA
  idxs.na.snps <- which( is.na(col_mean) )
  if (length(idxs.na.snps)>0) mat <- mat[, -idxs.na.snps] ##***
  
  dat.pca <- as.data.frame( mat )
  rownames(dat.pca) <- sampleID
  
  PCA <- prcomp(dat.pca)
  PC  <- predict(PCA)
  PC  <- data.frame(Sample_ID = rownames(PC), 
                    PC[, c("PC1", "PC2", "PC3")], 
                    stringsAsFactors = FALSE)
  
  write.table(PC, file = file.path(wk_dir, "result_PCA.txt"),
              quote = F, row.names = F, sep = "\t")
  
  ## plot PCA results
  p12 <- ggplot(data = PC, aes(PC1, PC2)) + 
    geom_point(size = 1) + 
    theme_bw() +
    theme(plot.title = element_text(size = 20, hjust = 0.5),
          axis.title = element_text(size = 15, face = "bold"),
          axis.text = element_text(size = 15, face = "bold")) + 
    ggtitle("PC2 ~ PC1")
  
  p13 <- ggplot(data = PC, aes(PC1, PC3)) + 
    geom_point(size = 1) + 
    theme_bw() +
    theme(plot.title = element_text(size = 20, hjust = 0.5),
          axis.title = element_text(size = 15, face = "bold"),
          axis.text = element_text(size = 15, face = "bold")) + 
    ggtitle("PC3 ~ PC1")
  
  p23 <- ggplot(data = PC, aes(PC2, PC3)) + 
    geom_point(size = 1) + 
    theme_bw() +
    theme(plot.title = element_text(size = 20, hjust = 0.5),
          axis.title = element_text(size = 15, face = "bold"),
          axis.text = element_text(size = 15, face = "bold")) + 
    ggtitle("PC3 ~ PC2")
  
  png(filename = file.path(wk_dir, "plots_PCA.png"),
      width = 12, height = 12, units = "in", res = 512)  
    p <- plot_grid(p12, p13, p23, nrow = 2)
    print(p)
  dev.off()


