#!/usr/bin/env Rscript

suppressMessages(library(optparse))

option_list = list(
  make_option(c("-t", "--type"), action = "store", type = "character", default = NA,
              help = "Job submission type (0 - initial submission, 1 - resubmission of failed jobs)"),
  make_option(c("-p", "--datapath"), action = "store", type = "character", default = NA,
              help = "Path to the directory containing necessary input data."),
  make_option(c("-o", "--resultpath"), action = "store", type = "character", default = NA,
              help = "Path to the directory for saving results."),
  make_option(c("-m", "--matrixpath"), action = "store", type = "character", default = NA,
              help = "Path to chromosome-wise LRR and BAF matrices."),
  make_option(c("-s", "--sourcefile"), action = "store", type = "character", default = NA,
              help = "Path to the scripts directory containing R scripts to be loaded into R."),
  make_option(c("-d", "--duplicates"), action = "store_true", default = FALSE,
              help = "[optional] Whether duplicate pairs information will be annotated in diagnosis plots."),
  make_option(c("-n", "--plot"), action = "store_true", default = FALSE,
              help = "[optional] Whether to generate diagnosis plots."),
  make_option(c("-r", "--script"), action = "store", type = "character", default = NA,
              help = "Path to the main script CNV.genotype.one.chr.one.batch.R."),
  make_option(c("-l", "--joblog"), action = "store", type = "character", default = NA,
              help = "path save all jobs log.")
)

opt = parse_args(OptionParser(option_list = option_list))
pars = c(opt$type, opt$datapath, opt$resultpath, opt$joblog,
         opt$matrixpath, opt$sourcefile, opt$script)

if ( any(is.na(pars)) ) {
  stop("All required parameters must be supplied. (--help for detail)")
}

script <- file.path(opt$script, "CNV.genotype.one.chr.one.batch.R")
cmd    <- paste("Rscript", script, 
                "--type", opt$type,
                "--datapath", opt$datapath,
                "--resultpath", opt$resultpath,
                "--matrixpath", opt$matrixpath,
                "--sourcefile", opt$sourcefile)

if ( opt$duplicates ) cmd <- paste(cmd, "--duplicates")
if ( opt$png ) cmd <- paste(cmd, "--plot")

path_joblog <- opt$joblog
if (!dir.exists(paths = path_joblog)) dir.create(path = path_joblog, showWarnings = F, recursive = T)
dir.create(path = file.path(path_joblog, "job", "ERROR"), showWarnings = F, recursive = T)
dir.create(path = file.path(path_joblog, "job", "OUT"), showWarnings = F, recursive = T)

path_job_error <- file.path(path_joblog, "job", "ERROR")
path_job_out   <- file.path(path_joblog, "job", "OUT")

file_cnvr <- "cnvr_batch.txt"  ## with batch information
dat_cnvr  <- read.delim(file = file.path(opt$datapath, file_cnvr), as.is = TRUE)
chrs = sort( unique(dat_cnvr$chr) )

for ( chr1 in chrs ) {
  
  dat_cnvr_chr1 = subset(dat_cnvr, chr == chr1)
  batch_chr1 = sort( unique(dat_cnvr_chr1$batch) )
  
  if ( nrow(dat_cnvr_chr1) == 0) {
    next
  }
  
  for ( batch1 in batch_chr1 ){
    
    cat("chr:", chr1, "batch1:", batch1, "\n")
    cmd1 = paste(cmd, "--chr", chr1, "--batch", batch1)
    bsub.cmd = paste("bsub -n 2 -W 10:00 -R 'rusage[mem=20000]' -P [account]", ## need to modify based on specific system
                     "-e", file.path(path_job_error, paste0("chr_", chr1, "_batch_", batch1, ".e")), 
                     "-o", file.path(path_job_out, paste0("chr_", chr1, "_batch_", batch1, ".o")),
                     "-q premium", shQuote(cmd1))
    cat(bsub.cmd, "\n")
    system(bsub.cmd)
    
    Sys.sleep(0.1)
  }
  
}



