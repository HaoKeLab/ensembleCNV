#!/usr/bin/env Rscript

suppressMessages(library(optparse))

option_list <- list(
  make_option(c("-o", "--resultpath"), action = "store", type = "character", default = NA,
              help = "path save all refinement results (cnvrs_refinement).")
)

opt <- parse_args(OptionParser(option_list = option_list))

pars <- c(opt$resultpath)
if (any(is.na(pars))) {
  stop("All parameters must be supplied. (--help for detail)")
}

path_result <- opt$resultpath
# combine refinement results ----------------------------------------------
path_refine <- file.path(path_result, "res_refine")
folders_chr1 <- list.files(path = path_refine, pattern = "^chr")

res_refinement <- data.frame()
for ( folder.chr1 in folders_chr1 ) {
  
  chr1 <- gsub("^chr", "", folder.chr1, perl = T)
  file.chr1 <- paste("CNVR_refine_chr_", chr1, "_detail.rds", sep = "")
  path.chr1.data <- file.path(path_refine, folder.chr1, "data")
  res.chr1 <- readRDS( file = file.path( path.chr1.data, file.chr1) )
  
  res_refinement <- rbind(res_refinement, res.chr1)
}

res_refinement$identicalID <- paste(res_refinement$Chr,
                                    res_refinement$snp.start.refine,
                                    res_refinement$snp.end.refine, sep = "___")

res_refinement_same <- subset(res_refinement, type.overlap.based.on.raw == "same")
res_refinement_refine <- subset(res_refinement, type.overlap.based.on.raw != "same")

res_refinement_refine <- subset(res_refinement_refine, 
                                !identicalID %in% res_refinement_same$identicalID)

# de-dulplicate CNVR
res_refinement_refine <- res_refinement_refine[!duplicated(res_refinement_refine$identicalID), ]
nrow(res_refinement_refine)

cnvrID_refine_same <- res_refinement_same$CNVR_ID

# clean -------------------------------------------------------------------

cnvrID_keep <- readRDS(file = file.path(path_result, "cnvrs_keep.rds"))

dat_cnvr_keep <- readRDS(file = file.path(path_result, "dat_cnvrs_keep.rds"))
dat_cnvr_keep$identicalID <- paste(dat_cnvr_keep$chr,
                                   dat_cnvr_keep$start_snp,
                                   dat_cnvr_keep$end_snp, sep = "___")

res_refinement_refine_clean <- subset( res_refinement_refine, !identicalID %in% dat_cnvr_keep$identicalID )

cnvrID_keep_final <- union(cnvrID_keep, cnvrID_refine_same)

saveRDS(cnvrID_final_keep, file = file.path(path_result, "cnvrs_keep_after_refine.rds"))
saveRDS(res_refinement_refine_clean, file = file.path(path_result, "dat_cnvrs_regenotype.rds"))


