


path_output <- ""

path_ensembleCNV <- ""
fileName_ensembleCNV <- ""

mat_ensembleCNV_before_refine <- readRDS( file = file.path(path_ensembleCNV, fileName_ensembleCNV) )
cnvrs   <- rownames(mat_ensembleCNV_before_refine)
samples <- colnames(mat_ensembleCNV_before_refine)

## keep 
path_res_refine <- ""
cnvrs_final_keep <- readRDS(file = file.path(path_res_refine, "cnvrs_final_keep.rds"))

mat_ensembleCNV_keep <- mat_ensembleCNV_before_refine[cnvrs_final_keep, ]

## refine
# after cnvrs need to regenotype, we can get matrix
path_res_regenotype <- ""
mat_ensembleCNV_regenotype <- readRDS(file = file.path(path_res_regenotype, "matrix_cnvrs_regenotype.rds"))

samples.regenotype <- colnames( mat_ensembleCNV_regenotype )

stopifnot( sum(samples.regenotype %in% samples) == length(samples))

mat_ensembleCNV_regenotype <- mat_ensembleCNV_regenotype[, samples]

## final
mat_ensembleCNV <- rbind(mat_ensembleCNV_keep, mat_ensembleCNV_regenotype)
dim(mat_ensembleCNV)

saveRDS(mat_ensembleCNV, file = file.path(path_output, "mat_ensembelCNV.rds"))









