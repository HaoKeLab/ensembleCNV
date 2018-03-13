
# input: matrix_cnvrs
plot_heatmap <- function(matrix_cnvrs, dt_cnvrs) {
  
  df_cnvrs <- as.data.frame(matrix_cnvrs)
  # cor
  cor1 <- cor(df_cnvrs, use = "na.or.complete")
  n <- ncol(cor1)
  # add columns annotate
  dt_cnvrs <- dt_cnvrs[order(dt_cnvrs$Sample_ID, dt_cnvrs$Position), ]
  dt_cnvr1 <- dt_cnvrs[1:n, ]
  snps_boundary <- ifelse(dt_cnvr1$snp_flag == 0, "add", "raw")
  snps_final    <- ifelse(dt_cnvr1$snp_flag == 2, "innerb", "outerb")
  
  annotate_col <- data.frame(
    group1 = snps_boundary,
    group2 = snps_final,
    stringsAsFactors = FALSE
  )
  rownames(annotate_col) <- dt_cnvr1$Name
  
  pheatmap(mat = cor1,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           annotation_col = annotation_col)
  
}