#!/usr/bin/env Rscript

suppressMessages(library(optparse))

option_list <- list(
  make_option(c("-p", "--datapath"), action = "store", type = "character", default = NA,
              help = "Path to the directory containing necessary input data."),
  make_option(c("-o", "--resultpath"), action = "store", type = "character", default = NA,
              help = "Path to the directory for saving results."),
  make_option(c("-m", "--matrixpath"), action = "store", type = "character", default = NA, 
              help = "Path to chromosome-wise LRR and BAF matrices."),
  make_option(c("-i", "--refinescript"), action = "store", type = "character", default = NA,
              help = "Path to the main script CNVR.boundary.refinement.R."),
  make_option(c("-s", "--rcppfile"), action = "store", type = "character", default = NA,
              help = "Path to refine.rcpp to be used in this R script."),
  make_option(c("-r", "--centromere"), action = "store", type = "character", default = NA,
              help = "Path to file with centromere position of each chromosome."),
  make_option(c("-n", "--plot"), action = "store_true", default = FALSE,
              help = "[optional] Whether to generate diagnosis plots.")
)

opt <- parse_args(OptionParser(option_list = option_list))
pars = c(opt$datapath, opt$resultpath, opt$matrixpath,
         opt$rcppfile, opt$centromere, opt$refinescript)

if ( any(is.na(pars)) ) {
  stop("All parameters must be supplied. (--help for detail)")
}

script_refine   <- opt$refinescript
path_result     <- opt$resultpath
path_matrix     <- opt$matrixpath
path_data       <- opt$datapath
script_rcpp     <- opt$rcppfile
file_centromere <- opt$centromere
flag_plot       <- opt$plot

# cnvrs refinement
dat_cnvrs_refine <- read.delim( file = file.path(path_result, "cnvr_refine.txt"), as.is = TRUE ) 
stopifnot( nrow(dat_cnvrs_refine) > 0 )

chrs <- sort( unique(dat_cnvrs_refine$chr))

cmd <- paste("Rscript", script_refine,
             "--datapath", path_data,
             "--resultpath", path_result,
             "--matrixpath", path_matrix,
             "--rcppfile", script_rcpp,
             "--centromere", file_centromere)

if ( flag_plot ) {
  cmd <- paste(cmd, "--plot")
} 

for (chr1 in chrs) {
  
  cmd.chr1 <- paste(cmd, 
                    "--chr", chr1)

  path_log <- file.path(path_result, "res_refine/chr", chr1, "log")
  if (!dir.exists(path_log)) dir.create(path = path_log, recursive = TRUE)

  bsub.cmd.chr1 <- paste("bsub -n 2 -W 10:00",
                         "-R 'rusage[mem=10000]'",
                         "-P [accout]",   ## need to modify based on specific system
                         "-J", paste0("chr", chr1),
                         "-q premium",
                         "-e", file.path(path_log, paste0("boundary_refine_chr", chr1, ".err")),
                         "-o", file.path(path_log, paste0("boundary_refine_chr", chr1, ".log")),
                         shQuote( cmd.chr1 ))
  
  cat("chr:", chr1, bsub.cmd.chr1, "\n")
  system( bsub.cmd.chr1 )
  Sys.sleep(0.1)
}
