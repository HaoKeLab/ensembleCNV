#!/usr/bin/env Rscript

## PCA
filename_matrix <- ""
path_input <- ""

suppressMessages({
  require(data.table)
  require(tibble)
  require(cowplot)
  require(ggplot2)
})

read_matrix <- function(path_input, filename_matrix) {
  
  dat_LRR = fread(input = file.path( path_input, filename_matrix) )
  
  dat_LRR <- as.data.frame(dat_LRR, stringsAsFactors = FALSE)
  dat_LRR <- column_to_rownames(dat_LRR, var = "V1")
  
  sampleID = rownames( dat_LRR )
  ## deal NA values in matrix
  mat = as.matrix(dat_LRR)
  rownames(mat) = sampleID
  colnames(mat) = NULL
  
  col_mean = colMeans(mat, na.rm = TRUE)
  for (i in 1:nrow(mat)) {
    v1 = as.vector(mat[i, ])
    idx1 = which(is.na(v1))
    if (length(idx1) >= 1) {
      mat[i, idx1] = col_mean[idx1]
    }
  }
  ## check which snps all values are NA
  idxs.na.snps <- which( is.na(col_mean) )
  mat.pca <- mat[, idxs.na.snps]
  
  dat.pca = as.data.frame( mat.pca )
  rownames(dat.pca) = sampleID
  
  PCA = prcomp(dat.pca)
  PC  = predict(PCA)
  PC  = as.data.frame(PC, stringsAsFactors = FALSE)
  PC  = PC[, c("PC1", "PC2", "PC3")]
  
  saveRDS(PC, file = file.path(path_input, "result_PCA.rds"))
  
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
  
  png(filename = "plots_PCA.png",
      width = 12, height = 12, units = "in", res = 512)
  p <- plot_grid(p12, p13, p23, nrow = 2)
  print(p)
  dev.off()
  
}


