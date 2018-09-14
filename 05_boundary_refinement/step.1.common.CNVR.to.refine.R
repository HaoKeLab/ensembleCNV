#!/usr/bin/env Rscript

suppressMessages(library(optparse))

option_list <- list(
  make_option(c("-p", "--datapath"), action = "store", type = "character", default = NA,
              help = "Path to the directory containing necessary input data."),
  make_option(c("-o", "--resultpath"), action = "store", type = "character", default = NA,
              help = "Path to the directory for saving results."),
  make_option(c("-c", "--freq"), action = "store", type = "double", default = NA,
              help = "Frequency cut-off to select CNVRs with common CNVs for boundary refinement.")
)


opt = parse_args(OptionParser(option_list = option_list))
pars = c(opt$datapath, opt$resultpath, opt$freq)

if ( any(is.na(pars)) ) {
  stop("All three parameters must be supplied. (--help for detail)")
}

cutoff_freq <- as.numeric( opt$freq )
path_data   <- opt$datapath
path_result <- opt$resultpath

path_output <- file.path( path_result ) ##"cnvr_refinement"
if (!dir.exists(paths = path_output) ) dir.create(path = path_output, showWarnings = F, recursive = T)

# the copy number matrix generated from CNV genotyping step
mat_CN <- readRDS( file = file.path(path_data, "matrix_CN.rds"))
n.sample <- ncol( mat_CN )
n.CNVR   <- nrow( mat_CN )

cnvrs <- rownames( mat_CN )

freqs_CNVR <- unlist( lapply(1:n.CNVR, FUN = function(i) {
  v1 <- as.integer( mat_CN[i, ])
  n1 <- sum( v1 %in% c(0, 1, 3))
  n1
}))

idxs.refine <- which( freqs_CNVR >= n.sample*cutoff_freq)

dat_freq <- data.frame(CNVR_ID = cnvrs, 
                       Freq = freqs_CNVR,
                       stringsAsFactors = F)

if (length(idxs.refine) > 0) {
  cnvrs_refine <- cnvrs[ idxs.refine ]
  cnvrs_keep   <- cnvrs[ -idxs.refine ]
} else {
  cnvrs_refine <- NULL
  cnvrs_keep   <- cnvrs
}

# write.table( data.frame(CNVR_ID = cnvrs_refine, stringsAsFactors = FALSE), 
             # file = file.path(path_output, "cnvrs_refine.txt"),
             # quote = F, row.names = F, col.names = F, sep = "\t")
# write.table( data.frame(CNVR_ID = cnvrs_keep, stringsAsFactors = FALSE),
             # file = file.path(path_output, "cnvrs_keep.txt"),
             # quote = F, row.names = F, col.names = F, sep = "\t")

file_cnvr <- "cnvr_genotype.txt"  ## with CNV genotype information
dat_cnvrs <- read.delim(file = file.path(path_data, file_cnvr), as.is = TRUE)
nms <- names(dat_cnvrs)
names(dat_cnvrs)[nms == "Freq"] <- "raw_Freq"
dat_cnvrs <- subset(dat_cnvrs, genotype == 1)

dat_cnvrs <- merge( dat_cnvrs, dat_freq, by = "CNVR_ID", all = FALSE)
dat_cnvrs <- dat_cnvrs[order(dat_cnvrs$chr, dat_cnvrs$posStart, dat_cnvrs$posEnd), ]
stopifnot( nrow(dat_cnvrs) == nrow(dat_freq) )

if (length(cnvrs_refine) > 0) {
  dat_cnvrs_refine <- subset( dat_cnvrs, CNVR_ID %in% cnvrs_refine )
  dat_cnvrs_keep   <- subset( dat_cnvrs, CNVR_ID %in% cnvrs_keep )
} else {
  dat_cnvrs_refine <- data.frame(NULL)
  dat_cnvrs_keep   <- dat_cnvrs
}

write.table( dat_cnvrs_keep, 
             file = file.path(path_output, "cnvr_keep.txt"),
             quote = F, row.names = F, sep = "\t")
write.table( dat_cnvrs_refine, 
             file = file.path(path_output, "cnvr_refine.txt"),
             quote = F, row.names = F, sep = "\t")



