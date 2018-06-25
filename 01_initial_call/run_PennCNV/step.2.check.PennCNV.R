#!/usr/bin/env Rscript

suppressMessages({
  require( optparse, quietly = TRUE)
})

options(warn = 2)

option_list <- list(
  make_option(c("-a", "--dat"), action = "store", default = NA, type = "character",
              help = "dat path from each sample for running PennCNV."),
  make_option(c("-b", "--main"), action = "store", default = NA, type = "character",
              help = "main path for each sample's list file and res."),
  make_option(c("-c", "--pfb"), action = "store", default = NA, type = "character",
              help = "pfb file."),
  make_option(c("-d", "--gcmodel"), action = "store", default = NA, type = "character",
              help = "gcmodel file."),
  make_option(c("-e", "--hmm"), action = "store", default = NA, type = "character",
              help = "hmm file.")
)

opt = parse_args(OptionParser(option_list = option_list))

path_dat     <- opt$dat
path_main    <- opt$main
file_pfb     <- opt$pfb
file_gcmodel <- opt$gcmodel
file_hmm     <- opt$hmm

if (any(is.na(c(path_dat, path_main, file_pfb, file_gcmodel, file_hmm)))) {
  stop("all parameters must be supplied. (--help for detail)")
}

path_list  <- file.path(path_main, "list")
path_res  <- file.path(path_main, "res")  ## PennCNV results folder

# submit jobs functions ---------------------------------------------------

cmd_PennCNV <- function(file_hmm, file_pfb, file_gcmodel, 
                        filename_sample, path_list, path_res_sample) {
  
  file_list <- file.path(path_list, filename_sample) 
  
  samplename <- gsub(pattern = ".txt$", replacement = "", filename_sample)
  file_log   <- file.path(path_res_sample, paste0(samplename, ".log"))
  file_rawcnv <- file.path(path_res_sample, paste0(samplename, ".rawcnv"))
  
  cmd <- paste("/hpc/packages/minerva-common/penncnv/2011Jun16/bin/detect_cnv.pl -test --confidence",
               "-hmm", file_hmm,
               "-pfb", file_pfb,
               "-gcmodel", file_gcmodel,
               "-list", file_list,
               "-log", file_log,
               "-out", file_rawcnv)
  
  cmd
} 

cmd_submitjob <- function(cmd.sample, samplename) {
  
  bsub.cmd <- paste("bsub -n 2 -W 00:30 -R 'rusage[mem=5000]' -P acc_haok01a", ##-R 'span[ptile=6]' 
                    "-J", samplename,
                    "-q premium",
                    shQuote(cmd))
  
  bsub.cmd
}

# main loop ---------------------------------------------------------------

sample_files <- list.files(path = path_dat)
cat("number of samples:", length(sample_files), "\n")

n.success <- 0
n.fail <- 0
for ( i in 1:length(sample_files) ) {
  
  sample_file <- sample_files[i]
  samplename  <- gsub(pattern = ".txt$", replacement = "", sample_file)
  
  path_res_sample <- file.path(path_res, samplename)
  file_rawcnv <- file.path(path_res_sample, paste0(samplename, ".rawcnv"))
  
  flag.folder <- dir.exists(paths = path_res_sample)
  flag.rawcnv <- file.exists(file_rawcnv)
  
  if ( flag.folder & flag.rawcnv ) {
    cat("Sample_ID:", samplename, "SUCCESS\n")
    n.success <- n.success + 1
  } else {
    
    cat("Sample_ID:", samplename, "FAILED\n")
    dir.create(path = path_res_sample, showWarnings = FALSE, recursive = TRUE)
    
    cmd.sample <- cmd_PennCNV(file_hmm = file_hmm,
                              file_pfb = file_pfb,
                              file_gcmodel = file_gcmodel,
                              filename_sample = sample_file,
                              path_list = path_list,
                              path_res_sample = path_res_sample)
    
    cmd.job   <- cmd_submitjob(cmd.sample = cmd.sample, samplename = samplename)
    
    system(cmd.job)
    Sys.sleep(0.1)
    
    n.fail <- n.fail + 1
    
  }
}

cat("summary, total:", length(sample_files),
    "number of success:", n.success,
    "number of fail:", n.fail, "\n")






