#!/usr/bin/env Rscirpt --vanilla

args <- commandArgs(trailingOnly = TRUE)

file_cnv    <- args[1]
cnvCaller   <- args[2]
path_output <- args[3]

generate_matrix_one_method <- function(dt_cnvs) {
  
  dt_cnvs$identical_ID <- paste(dt_cnvs$Sample_ID,
                                dt_cnvs$CNVR_ID,
                                sep = "_")
  
  samples <- unique(dt_cnvs$Sample_ID)
  cnvrs <- unique(dt_cnvs$CNVR_ID)
  
  n.samples <- length(samples)
  n.cnvrs <- length(cnvrs)
  
  mat <- matrix(rep(NA, n.samples*n.cnvrs),
                nrow = n.cnvrs, ncol = n.samples)
  
  rownames(mat) <- cnvrs
  colnames(mat) <- samples
  
  for ( i in 1:n.cnvrs ) {
    
    cnvr1 <- cnvrs[i]
    cat("idx:", i, "in", n.cnvrs, "CNVR_ID:", cnvr1, "\n")
    
    dt_cnvs1 <- subset(dt_cnvs, CNVR_ID == cnvr1)  ## all 
    tbl1 <- table(dt_cnvs1$identical_ID)
    
    idx1 <- which(tbl1 == 2)
    if (length(idx1) >= 1) {
      id1 <- names(tbl1)[idx1]
      dt_cnvs1 <- subset(dt_cnvs1, !identical_ID %in% id1) ## exclude confict CNVR_ID_Sample_ID
    }
    
    idxs.match <- match(samples, dt_cnvs1$Sample_ID)
    mat[i, ] <- dt_cnvs1$CN[idxs.match]
    
  }
  
  mat[which(is.na(mat))] <- 2
  return(mat)
}

dt_cnvs <- readRDS(file = file_cnv)
names(dt_cnvs)[which(names(dt_cnvs) == "CNV_Value")] <- "CN"
mat <- generate_matrix_one_method(dt_cnvs = dt_cnvs)

saveRDS(mat, file = file.path(path_output, paste0("matrix_", cnvCaller, ".rds")))






