#!/hpc/packages/minerva-common/R/3.3.1/lib64/R/bin/Rscript --vanilla

## all data save as vector to fast speed

# create folder name with clean CNVR data
file_cnvr <- "cnvrs_batch_annotate.rds"  # with batch annotated
path_cnvr <- "/sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/ensembleCNV/5batch/res"

dt_cnvr <- readRDS(file = file.path(path_cnvr, file_cnvr))

tbl <- table(dt_cnvr$chr, dt_cnvr$batch)

dat <- as.data.frame(tbl)
names(dat) <- c("chr", "batch", "Freq") ##

dat1 <- subset(dat, Freq != 0)
cat("total number of folders:", nrow(dat1), "\n") #
# before doing this please check if all chr and batch finished
# path_res_chr_batch

fun_combine <- function(path_res, n_sample, path_cnvr, file_cnvr) {
  
  dt_cnvr <- readRDS(file = file.path(path_cnvr, file_cnvr))
  tbl <- table(dt_cnvr$chr, dt_cnvr$batch)
  
  n_cnvr <- sum(tbl)
  dat <- as.data.frame(tbl)
  names(dat) <- c("chr", "batch", "Freq")
  
  dat <- subset(dat, Freq != 0)
  
  cnvrs <- character(length = n_cnvr)
  samples <- character(length = n_sample)
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
      samples <- dat_cnvr1$Sample_ID

      if(all(copy_numbers == 2)) {
        idxs_del[nth_cnvr] <- 1
      }
      
      nth_cnvr <- nth_cnvr + 1
    }
    
  }
  
  return(list(CNVR_ID = cnvrs,
              Sample_ID = samples,
              Copy_numbers = copy_numbers,
              GQs = GQs,
              idxs_del = idxs_del))
}

# main ----------------------------------------------------------------
path_res_chr_batch <- "/sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/ensembleCNV/5batch/res_regenotypeCNV_raw/pred"

list_res <- fun_combine(path_res = path_res_chr_batch,
                        path_cnvr = path_cnvr, 
                        n_sample = 2765,  ## 
                        file_cnvr = file_cnvr)

path_output <- "/sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/ensembleCNV/5batch/res_regenotypeCNV_raw/res_combined"
saveRDS(list_res$CNVR_ID, file = file.path(path_output, "CNVR_ID.rds"))
saveRDS(list_res$Sample_ID, file = file.path(path_output, "Sample_ID.rds"))
saveRDS(list_res$Copy_numbers, file = file.path(path_output, "Copy_Number.rds"))
saveRDS(list_res$GQs, file = file.path(path_output, "GQ_score.rds"))
saveRDS(list_res$idxs_del, file = file.path(path_output, "idxs_del.rds"))


