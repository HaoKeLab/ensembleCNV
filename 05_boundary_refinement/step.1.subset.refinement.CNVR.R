

## select CNVR using cutoff_freq

path_ensembleCNV <- ""
fileName_ensembleCNV <- ""
fileName_CNVR <- ""
cutoff_freq <- 0.01
path_output <- ""

mat_ensembleCNV <- readRDS( file = file.path(path_ensembleCNV, fileName_ensembleCNV) )
n.sample <- ncol( mat_ensembleCNV )
n.CNVR   <- nrow( mat_ensembleCNV )

cnvrIDs  <- rownames( mat_ensembleCNV )

## calculate freq of CNVR

freq_CNVR <- unlist( lapply(1:nrow(mat_ensembleCNV), FUN = function(i) {
  v1 <- as.integer( mat_ensembleCNV[i, ] )
  n1 <- sum( v1 %in% c(0, 1, 3) )
}))

idxs.refine <- which( freq_CNVR >= n.sample*cutoff_freq )
length(idxs.refine)

cnvrs_refine <- cnvrIDs[ idxs.refine ]
cnvrs_keep   <- cnvrIDs[ -idxs.refine ]

saveRDS( cnvrs_refine, file = file.path( path_output, "cnvrs_refine.rds") )
saveRDS( cnvrs_keep, file = file.path( path_output, "cnvrs_keep.rds"))

dat_cnvr <- readRDS(file = file.path(path_ensembleCNV, fileName_CNVR))

dat_cnvr_keep <- subset( dat_cnvr, CNVR_ID %in% cnvrs_keep )
dat_cnvr_refine <- subset( dat_cnvr, CNVR_ID %in% cnvrs_refine )

saveRDS( dat_cnvr_keep, file = file.path(path_output, "dat_cnvrs_keep.rds") )
saveRDS( dat_cnvr_refine, file = file.path(path_output, "dat_cnvrs_refine.rds"))










