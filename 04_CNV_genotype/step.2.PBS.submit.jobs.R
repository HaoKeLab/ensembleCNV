#!/usr/bin/env Rscript

## NOTE: The scripts embraced by "##<<<... ##>>>..." need to be specified based on your system

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
              help = "Path to the directory saving job logs.")
)

opt = parse_args(OptionParser(option_list = option_list))
pars = c(opt$type, opt$datapath, opt$resultpath, opt$joblog,
         opt$matrixpath, opt$sourcefile, opt$script)

if ( any(is.na(pars)) ) {
  stop("All required parameters must be supplied. (--help for detail)")
}

# script <- file.path(opt$script, "CNV.genotype.one.chr.one.batch.R")

# modified 
cmd <- paste0(  "-v TYPE=", opt$type,
                ",DTPATH=", opt$datapath,
                ",RESPATH=", opt$resultpath,
                ",MATPATH=", opt$matrixpath,
                ",SOURCE=", opt$sourcefile, 
                " ", file.path(opt$script, "run.CNV.genotype.sh") )


if ( opt$duplicates ) cmd <- paste(cmd, "--duplicates")
if ( opt$plot ) cmd <- paste(cmd, "--plot")

path_joblog <- opt$joblog
if (!dir.exists(paths = path_joblog)) dir.create(path = path_joblog, showWarnings = F, recursive = T)
dir.create(path = file.path(path_joblog, "job", "ERROR"), showWarnings = F, recursive = T)
dir.create(path = file.path(path_joblog, "job", "OUT"), showWarnings = F, recursive = T)

path_job_error <- file.path(path_joblog, "job", "ERROR")
path_job_out   <- file.path(path_joblog, "job", "OUT")

file_cnvr <- "cnvr_batch.txt"  ## with batch information
dat_cnvr  <- read.delim(file = file.path(opt$datapath, file_cnvr), as.is = TRUE)
chrs <- sort( unique(dat_cnvr$chr) )

for ( chr1 in chrs ) {
  
  dat_cnvr_chr1 = subset(dat_cnvr, chr == chr1)
  batch_chr1 = sort( unique(dat_cnvr_chr1$batch) )
  
  if ( nrow(dat_cnvr_chr1) == 0) {
    next
  }
  
  for ( batch1 in batch_chr1 ){
    
    cat("chr:", chr1, "batch1:", batch1, "\n")
    cmd1 = paste(cmd, "--chr", chr1, "--batch", batch1)

    qsub.log.file <- file.path(path_job_error, paste0("chr_", chr1, "_batch_", batch1, ".e"))
    qsub.err.file <- file.path(path_job_out, paste0("chr_", chr1, "_batch_", batch1, ".o"))
    
##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
## configure based on your system
    
    qsub.cmd <- paste("qsub -l nodes=2,ncpus=10,mem=10gb,walltime=0:15:00 -q <queque.name> -o", 
    qsub.log.file,
    "-e", qsub.err.file , 
    cmd1)

#    bsub.cmd = paste("bsub -n 2 -W 10:00 -R 'rusage[mem=20000]' -P <account>",
#                     "-e", file.path(path_job_error, paste0("chr_", chr1, "_batch_", batch1, ".e")), 
#                     "-o", file.path(path_job_out, paste0("chr_", chr1, "_batch_", batch1, ".o")),
#                     "-q premium", shQuote(cmd1))
##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    cat(qsub.cmd, "\n")
    system(qsub.cmd)
    
    Sys.sleep(0.1)
  }
  
}



