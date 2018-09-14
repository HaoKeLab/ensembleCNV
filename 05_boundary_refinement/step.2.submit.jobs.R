#!/usr/bin/env Rscript

suppressMessages(library(optparse))

option_list <- list(
  make_option(c("-i", "--refinescript"), action = "store", type = "character", default = NA,
              help = "script to refine CNVR."),
  make_option(c("-m", "--matrixpath"), action = "store", type = "character", default = NA, 
              help = "path save all LRR and BAF chromosome based data"),
  make_option(c("-n", "--plot"), action = "store_true", default = FALSE,
              help = "plot png or not"),
  make_option(c("-p", "--datapath"), action = "store", type = "character", default = NA,
              help = "path contain all running needed data"),
  make_option(c("-o", "--resultpath"), action = "store", type = "character", default = NA,
              help = "path save all results"),
  make_option(c("-s", "--rcppfile"), action = "store", type = "character", default = NA,
              help = "rcpp script need to be soucred"),
  make_option(c("-r", "--centromere"), action = "store", type = "character", default = NA,
              help = "centromere position of each chromosome.")
)

opt <- parse_args(OptionParser(option_list = option_list))
pars = c(opt$datapath, opt$resultpath, opt$matrixpath,
         opt$rcppfile, opt$centromere, opt$refinescript)

if ( any(is.na(pars)) ) {
  stop("All parameters must be supplied. (--help for detail)")
}

script_refine <- opt$refinescript
path_result <- opt$resultpath
path_matrix <- opt$matrixpath
path_data   <- opt$datapath
script_rcpp <- opt$rcppfile
file_centromere <- opt$centromere
flag_plot   <- opt$plot

dat_cnvrs_refine <- readRDS(file = file.path(path_result, "dat_cnvrs_refine.rds"))
chrs <- sort( unique(dat_cnvrs_refine$chr))

cmd <- paste(script_refine,
             "--matrixpath", path_matrix,
             "--datapath", path_data,
             "--resultpath", path_result,
             "--rcppfile", script_rcpp,
             "--centromere", file_centromere)
if ( flag_plot ) {
  cmd <- paste(cmd, "--plot")
} 


for (chr1 in chrs) {
  
  cmd.chr1 <- paste(cmd, 
                    "--chr", chr1)
  bsub.cmd.chr1 <- paste("bsub -n 2 -W 10:00",
                         "-R 'rusage[mem=10000]'",
                         "-P [accout]",
                         "-J", paste0("chr", chr1),
                         "-q premium",
                         shQuote( cmd.chr1 ))
  
  cat("chr:", chr1, bsub.cmd.chr1, "\n")
  system( bsub.cmd.chr1 )
  Sys.sleep(0.5)
}
