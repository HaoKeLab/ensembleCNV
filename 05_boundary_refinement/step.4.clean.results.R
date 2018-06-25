#!/usr/bin/env Rscript

path_refinement_result <- ""
path_output <- ""
# combine refinement results ----------------------------------------------

res_refinement <- data.frame()

for ( chr1 in 1:22 ) {
  
  cat("chr:", chr1, "\n") 
  
  folder.chr1 <- paste("chr", chr1, sep = "")
  if (!dir.exists(paths = file.path(path_refinement_result, folder.chr1))) next
  file.chr1   <- paste("CNVR_refine_chr_", chr1, "_detail.rds", sep = "")
  
  res.chr1 <- readRDS( file = file.path( path_refinement_result, folder.chr1, file.chr1) )
  
  res_refinement <- rbind(res_refinement, res.chr1)
  # type.overlap.based.on.raw
}

res_refinement$identicalID <- paste(res_refinement$Chr,
                                    res_refinement$snp.start.refine,
                                    res_refinement$snp.end.refine, sep = "___")

res_refinement_same <- subset(res_refinement, type.overlap.based.on.raw == "same")
res_refinement_refine <- subset(res_refinement, type.overlap.based.on.raw != "same")

res_refinement_refine <- subset(res_refinement_refine, 
                                !identicalID %in% res_refinement_same$identicalID)

## de-dulplicate CNVR
res_refinement_refine <- res_refinement_refine[!duplicated(res_refinement_refine$identicalID), ]
nrow(res_refinement_refine)

cnvrID_refine_same <- res_refinement_same$CNVR_ID

# clean -------------------------------------------------------------------

cnvrID_keep <- readRDS(file = file.path(path_output, "cnvrs_keep.rds"))

dat_cnvr_keep <- readRDS(file = file.path(path_output, "dat_cnvrs_keep.rds"))
dat_cnvr_keep$identicalID <- paste(dat_cnvr_keep$chr,
                                 dat_cnvr_keep$start_snp,
                                 dat_cnvr_keep$end_snp, sep = "___")

res_refinement_refine_clean <- subset(res_refinement_refine, !identicalID %in% dat_cnvr_keep$identicalID)
nrow(res_refinement_refine_clean)

cnvrID_final_keep <- union(cnvrID_keep, cnvrID_refine_same)

saveRDS(cnvrID_final_keep, file = file.path(path_output, "cnvrs_final_keep.rds"))
saveRDS(res_refinement_refine_clean, file = file.path(path_output, "dat_cnvrs_regenotype.rds"))



