
# plot heatmap

plot_heatmap <- function(dt_lrr_heatmap, dt_snps_flag) {
  
  cor1 <- cor(dt_lrr_heatmap, use = "na.or.complete")
  
  groups1 <- ifelse(dt_snps_flag$snp_flag == 0, "snps_add", "snps_raw")
  groups2 <- ifelse(dt_snps_flag$snp_flag == 2, "inner_boundary", "outer_boundary")
  annotation_col1 <- data.frame(
    group1 = groups1,
    group2 = groups2
  )
  rownames(annotation_col1) <- colnames(dt_lrr_heatmap)
  
  pheatmap(cor1, 
           cluster_cols = FALSE,
           cluster_rows = FALSE,
           annotation_col = annotation_col1)
  
}
