

path_output = ""
filename_output = ""
## summary CNVR_ID for each method
path_mat = "/sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/ensembleCNV.IPQ.each/5batch/iPattern"
filename_mat = "matrix_CNVR_sampleID_iPattern.rds"

path_mat = "/sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/ensembleCNV.IPQ.each/5batch/PennCNV"
filename_mat = "matrix_CNVR_sampleID_PennCNV.rds"
# 012 0123   02  023   12  123   23
# 987  687   72  124 7411 2446 2759

path_mat = "/sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/ensembleCNV.IPQ.each/5batch/QuantiSNP"
filename_mat = "matrix_CNVR_sampleID_QuantiSNP.rds"
# 012  0123    02   023    12   123    23
# 715  1635   126   133 10634  5423  6419

path_mat = "/sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/ensembleCNV/5batch/res"
filename_mat  = "matrix_ensembleCNV_after_regenotypeCNV.rds"
# 0    01   012  0123   013    02   023    03     1    12   123    13    23
# 11    62  1767   476    57   124    14     7    21 11044  3306    56  3438
# 3
# 24

mat = readRDS( file = file.path(path_mat, filename_mat) )
dim(mat)

cnvrID = rownames(mat)

dat_summary = data.frame()
for ( i in 1:nrow(mat) ) {
  
  cnvr1 = cnvrID[i]
  
  v1 = as.vector(mat[cnvr1, ])
  
  idx.nocall = which(v1 == -9)
  if (length(idx.nocall) >= 1) {
    v1 = v1[-idx.nocall]
  }
  
  cn_uniq = sort(unique(v1))
  cn_flag = paste(cn_uniq, collapse = "")
  
  v1f = factor(v1, levels = c(0 , 1, 2, 3))
  tbl1 = table(v1f)
  
  dat_cnvr1 = data.frame(CNVR_ID = cnvr1, flag = cn_flag,
                         n_CN_0 = tbl1[["0"]], n_CN_1 = tbl1[["1"]],
                         n_CN_2 = tbl1[["2"]], n_CN_3 = tbl1[["3"]],
                         stringsAsFactors = FALSE)
  
  dat_summary = rbind(dat_summary, dat_cnvr1)
}

table(dat_summary$flag)

write.table(dat_summary, file = file.path(path_output, filename_output),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)






