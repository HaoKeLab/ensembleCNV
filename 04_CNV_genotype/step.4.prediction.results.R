#!/usr/bin/env Rscript --vanilla

args <- commandArgs(trailingOnly = TRUE)

n <- as.integer( args[1] ) ## number of samples
path_cnvr <- args[2]
path_pred <- args[3]
path_res  <- args[4]

## all data save as vector to fast speed
n.samples <- n

# create folder name with clean CNVR data
file_cnvr <- "cnvrs_annotated_batch.rds"  # with batch annotated
dt_cnvr <- readRDS(file = file.path(path_cnvr, file_cnvr))

tbl <- table(dt_cnvr$chr, dt_cnvr$batch)

dat <- as.data.frame(tbl)
names(dat) <- c("chr", "batch", "Freq") ##

dat1 <- subset(dat, Freq != 0)
cat("total number of folders:", nrow(dat1), "\n") #

## cnvrs chr-batch-based results
path_res_chr_batch <- path_pred
path_output <- path_res

dt_cnvr <- readRDS(file = file.path(path_cnvr, file_cnvr))

tbl <- table(dt_cnvr$chr, dt_cnvr$batch)

dat <- as.data.frame(tbl)
names(dat) <- c("chr", "batch", "Freq") ##

dat1 <- subset(dat, Freq != 0)
cat("total number of folders:", nrow(dat1), "\n") #

fun_combine <- function(path_res, n_sample, path_cnvr, file_cnvr) {
  
  dt_cnvr <- readRDS(file = file.path(path_cnvr, file_cnvr))
  tbl <- table(dt_cnvr$chr, dt_cnvr$batch)
  
  n_cnvr <- sum(tbl)
  dat <- as.data.frame(tbl)
  names(dat) <- c("chr", "batch", "Freq")
  
  dat <- subset(dat, Freq != 0)
  
  cnvrs <- character(length = n_cnvr)
  samples <- character(length = n_sample)
  flag.samples <- FALSE
  copy_numbers <- integer(length = (n_cnvr*n_sample))
  GQs <- numeric(length = (n_cnvr*n_sample))
  idxs_del <- rep(0, n_cnvr)  # del CNVR raw(all CN = 2)
  
  nth_cnvr <- 1
  for (i in 1:nrow(dat)) {
    
    chr1 <- dat$chr[i]
    batch1 <- dat$batch[i]
    
    foldername1 <- paste0("chr_", chr1, "_batch_", batch1)
    path1 <- file.path(path_res, foldername1) # sub path
    
    cat("combine CNVR in the folder:", foldername1, "\n")
    
    files <- list.files(path = path1)
    for (k in 1:length(files)) {
      
      file1 <- files[k]
      dat_cnvr1 <- readRDS(file = file.path(path1, file1))
      dat_cnvr1 <- dat_cnvr1[order(dat_cnvr1$Sample_ID), ]
      
      cnvr1 <- unique(dat_cnvr1$CNVR_ID)
      
      idx <- ((nth_cnvr-1)*n_sample+1):(nth_cnvr*n_sample)
      copy_numbers[idx] <- dat_cnvr1$CN_gatk_pred ## output
      GQs[idx] <- dat_cnvr1$value_GQ  ##
      
      cnvrs[nth_cnvr] <- cnvr1
      # samples <- dat_cnvr1$Sample_ID
      if ( !flag.samples ) {
        samples <- dat_cnvr1$Sample_ID
      } 
      stopifnot( all(samples == dat_cnvr1$Sample_ID))
      
      if(all(copy_numbers == 2)) {
        idxs_del[nth_cnvr] <- 1
      }
      
      nth_cnvr <- nth_cnvr + 1
    }
    
  }
  
  mat_CN <- matrix(copy_numbers, nrow = length(cnvrs),
                   ncol = length(samples), dimnames = list(cnvrs, samples))
  mat_GQ <- matrix(GQs, nrow = length(cnvrs),
                   ncol = length(samples), dimnames = list(cnvrs, samples))
  
  mat_CN_clean <- mat_CN[!idxs_del == 1, ]
  mat_GQ_clean <- mat_GQ[!idxs_del == 1, ]
  
  cnvrs_clean  <- cnvrs[!idxs_del == 1]
  
  return(list(CNVR_ID = cnvrs,
              Sample_ID = samples,
              mat_CN = mat_CN_clean,
              mat_GQ = mat_GQ_clean))
  
}


list_res <- fun_combine(path_res = path_res_chr_batch,
                        path_cnvr = path_cnvr, 
                        n_sample = n.samples,  ## 
                        file_cnvr = file_cnvr)

saveRDS(list_res$CNVR_ID, file = file.path(path_output, "CNVR_ID.rds"))
saveRDS(list_res$Sample_ID, file = file.path(path_output, "Sample_ID.rds"))
saveRDS(list_res$mat_CN, file = file.path(path_output, "mat_CN.rds"))
saveRDS(list_res$mat_GQ, file = file.path(path_output, "mat_GQ.rds"))




