

# CNVR keep ---------------------------------------------------------------

## all CNVR after regenotype
path_ensembleCNV = "/sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/ensembleCNV/5batch/res"
file_ensembleCNV = file.path(path_ensembleCNV, "matrix_ensembleCNV_after_regenotypeCNV.rds")
mat_ensembleCNV = readRDS( file = file_ensembleCNV )
dim(mat_ensembleCNV)
n_sample = ncol(mat_ensembleCNV)
n_CNVR = nrow(mat_ensembleCNV)
CNVRs = rownames(mat_ensembleCNV)

## calculate set cutoff of freq (1%)
freq_CNVR = unlist( lapply(1:nrow(mat_ensembleCNV), FUN = function(i) {
  v1 = as.integer(mat_ensembleCNV[i, ])
  n1 = sum(v1 %in% c(0, 1, 3))
}) )

length(freq_CNVR)

index = which(freq_CNVR >= n_sample*0.01)
length(index) ## 3886 (number of CNVR)

cnvrs_refine = CNVRs[index]  ## refine
cnvrs_keep = CNVRs[-index]   ## keep

path_output = "/sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/ensembleCNV/after_boundary_refine/dat"
saveRDS(cnvrs_keep, file = file.path(path_output, "cnvrs_keep.rds"))
saveRDS(cnvrs_refine, file = file.path(path_output, "cnvrs_refine_raw.rds"))

length(cnvrs_keep)   ## 16521
length(cnvrs_refine) ## 3886



# ## submit jobs ----------------------------------------------------------

script = "/sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/boundary_refinement/cpp_version/boundary_refinement_V2_rcpp.R"

for (chr1 in 1:22) {
  
  cmd.chr1 = paste(script, chr1)
  bsub.cmd <- paste("bsub -n 2 -W 10:00 -R 'rusage[mem=10000]' -P acc_haok01a", ##-R 'span[ptile=6]'
                    "-J", chr1,
                    "-q premium",
                    shQuote(cmd.chr1))
  
  cat(bsub.cmd, "\n")
  system(bsub.cmd)
  Sys.sleep(0.5)
}

## 
# cd /sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/boundary_refinement/cpp_version
# module load R
# ./boundary_refinement_V2_rcpp.R 6(chr1)



# merge all results -------------------------------------------------------

path_res = "/sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/boundary_refinement/res_V2_rcpp"

res = data.frame()
for (chr1 in 1:22) {
  
  cat("chr:", chr1, "\n")
  
  folder.chr1 = paste0("chr", chr1)
  file.chr1 = paste0("CNVR_refine_chr_", chr1, "_detail.rds")
  
  res.chr1 = readRDS(file = file.path(path_res, folder.chr1, file.chr1))
  ## add flag
  res.chr1$identicalID = paste(res.chr1$Chr, 
                               res.chr1$snp.start.refine,
                               res.chr1$snp.end.refine, sep = "___")
  res.chr1$flag.identicalID = 0
  tbl.chr1 = table(res.chr1$identicalID)
  
  idx.chr1.dup = which(tbl.chr1 >= 2)
  if ( length(idx.chr1.dup) >= 1 ) {
    
    nms.chr1 = names(tbl.chr1)[idx.chr1.dup]
    res.chr1$flag.identicalID[which(res.chr1$identicalID %in% nms.chr1)] = 1
  }
  
  res = rbind(res, res.chr1)  ## final refine result
}

table(res$type.overlap.based.on.raw)
nrow(res)

stopifnot(nrow(res) == 3886)

# clean res ---------------------------------------------------------------

nrow(res)

identicalID.unique = unique(res$identicalID)
res_clean = data.frame()
for ( i in 1:length(identicalID.unique) ) {
  
  id1 = identicalID.unique[i]
  dat1= subset(res, identicalID == id1)
  
  if (nrow(dat1) == 1) {
    res_clean = rbind(res_clean, dat1)
  } else {
    
    type_dat1 = dat1$type.overlap.based.on.raw
    idx1 = which(type_dat1 == "same") 
    
    if (length(idx1) >= 1) {
      dat11 = dat1[idx1[1], ]
      res_clean = rbind(res_clean, dat11)
    } else {
      dat11 = dat1[1, ]
      res_clean = rbind(res_clean, dat11)
    }
    
  }
  
}

nrow(res)
nrow(res_clean)
tbl_clean = table(res_clean$identicalID)
stopifnot(all(tbl_clean == 1))

path_output = "/sc/orga/projects/haok01a/chengh04/paper/ensembleCNV/FA/data/res_regenotype"
res_del = subset(res, !CNVR_ID %in% res_clean$CNVR_ID)
nrow(res_del)


saveRDS(res_del$CNVR_ID, file = file.path(path_output, "cnvrs_refine_del.rds"))

table(res_clean$type.overlap.based.on.raw)

names(res_clean)[which(names(res_clean) == "Chr")] = "chr"
names(res_clean)[which(names(res_clean) == "snp.start.refine")] = "start_snp"
names(res_clean)[which(names(res_clean) == "snp.end.refine")] = "end_snp"

res_clean_refine = subset(res_clean, type.overlap.based.on.raw != "same")
nrow(res_clean_refine)

res_clean_same =  subset(res_clean, type.overlap.based.on.raw == "same")
nrow(res_clean_same)


for ( i in 1:nrow(res_clean_refine) ) {
  
  id1 = res_clean_refine$identicalID[i]
  idx1 = which(res_clean_same$identicalID == id1)
  
  cat(id1, ifelse(length(idx1) == 1, "delete", "keep"), "\n")
}


# path_output = "/sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/boundary_refinement/res_V2_rcpp"
file_cnvr_clean = "cnvrs_refine.rds"
path_cnvr_clean = "/sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/ensembleCNV/after_boundary_refine/dat"

saveRDS(res_clean_refine, file = file.path(path_cnvr_clean, "cnvrs_refine.rds"))
saveRDS(res_clean_same, file = file.path(path_cnvr_clean, "cnvrs_refine_same.rds"))


# generate same cnvr ------------------------------------------------------

path_cnvr_clean = "/sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/ensembleCNV/after_boundary_refine/dat"

## get these CNVR results from first time gmm modeling
cnvrs_keep = readRDS(file = file.path(path_cnvr_clean, "cnvrs_keep.rds"))
cnvrs_refine_same = readRDS(file = file.path(path_cnvr_clean, "cnvrs_refine_same.rds"))
cnvrs_refine_same = cnvrs_refine_same$CNVR_ID


cnvrs_keep_all = c(cnvrs_keep, cnvrs_refine_same) ## 
# simple way direct from ensembleCNV matrix
path_ensembleCNV = "/sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/ensembleCNV/5batch/res"
file_ensembleCNV = file.path(path_ensembleCNV, "matrix_ensembleCNV_after_regenotypeCNV.rds")
mat_ensembleCNV = readRDS( file = file_ensembleCNV )

## keep_all part mat_ensembleCNV
mat_ensembleCNV_keep = mat_ensembleCNV[cnvrs_keep_all, ]
dim(mat_ensembleCNV_keep)

saveRDS(mat_ensembleCNV_keep, file = file.path(path_cnvr_clean, "mat_ensembleCNV_keep_refine.rds"))
































