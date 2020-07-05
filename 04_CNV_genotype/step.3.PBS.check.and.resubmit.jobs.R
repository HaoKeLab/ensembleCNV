#!/usr/bin/env Rscript

## NOTE: The scripts embraced by "##<<<... ##>>>..." need to be specified based on your system

suppressMessages(require(optparse))

option_list = list(
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
              help = "Path to the directory saving job logs."),
  make_option(c("-f", "--flag"), action = "store", type = "integer", default = NA,
              help = "0: only print the running status of CNV genotyping; 1: resubmit jobs for unfinished CNV genotyping")
)

opt = parse_args(OptionParser(option_list = option_list))
pars = c(opt$type, opt$datapath, opt$resultpath, opt$joblog,
         opt$matrixpath, opt$sourcefile, opt$script, opt$flag)

if ( any(is.na(pars)) ) {
  stop("All parameters must be supplied. (--help for detail)")
}

flag <- as.integer( opt$flag )  ## 0 or 1

# resubmit unfinished jobs
file_cnvr <- "cnvr_batch.txt"  ## with batch information
dt_cnvr_raw <- read.delim(file = file.path(opt$datapath, file_cnvr), as.is = TRUE)
dt_cnvr_raw <- dt_cnvr_raw[order(dt_cnvr_raw$chr, dt_cnvr_raw$batch), ]
# add fname column
dt_cnvr_raw$fname <- paste0(dt_cnvr_raw$CNVR_ID, "_pred.rds") 

tbl_raw <- table(dt_cnvr_raw$chr, dt_cnvr_raw$batch)
dt_freq_raw <- as.data.frame(tbl_raw)
names(dt_freq_raw) <- c("chr", "batch", "Freq")

dt_freq_raw <- subset(dt_freq_raw, Freq != 0)  ## subset non-null batch
dt_freq_raw <- dt_freq_raw[order(dt_freq_raw$chr, dt_freq_raw$batch), ]

path_main_pred <- file.path(opt$resultpath, "pred")
path_main_failed <- file.path(opt$resultpath, "cnvrs_error")

# create script
script <- file.path(opt$script, "CNV.genotype.one.chr.one.batch.R")

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

# check if CNV genotyping for all CNVRs is finished ----------------------------------
check_jobs <- function(path_main, dt_cnvr_raw, flag, path_main_failed, path_joblog) {
  
  path_job_error <- file.path(path_joblog, "job", "ERROR")
  path_job_out   <- file.path(path_joblog, "job", "OUT")
  
  # remove all previous results
  system( paste("rm -rf", path_main_failed) )
  
  tbl_raw <- table(dt_cnvr_raw$chr, dt_cnvr_raw$batch)
  dt_freq_raw <- as.data.frame(tbl_raw)
  names(dt_freq_raw) <- c("chr", "batch", "Freq")
  
  dt_freq_raw <- subset(dt_freq_raw, Freq != 0)
  dt_freq_raw <- dt_freq_raw[order(dt_freq_raw$chr, dt_freq_raw$batch), ]
  
  for (i in 1:nrow(dt_freq_raw)) {
    
    chr1 <- dt_freq_raw$chr[i]
    batch1 <- dt_freq_raw$batch[i]
    freq1 <- dt_freq_raw$Freq[i]
    
    foldername1 <- paste0("chr_", chr1, "_batch_", batch1)
    path1 <- file.path(path_main, foldername1)
    
    if ( !dir.exists(paths = path1) ) { 
      cat("CHR:", chr1, "BATCH:", batch1, "The whole batch failed and jobs will be resubmitted.\n")
      
      # submit jobs
      if (flag == 1) {
        
        cmd1 = paste(cmd, "--chr", chr1, "--batch", batch1, "--type", 0)

##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
## configure based on your system

    qsub.cmd <- paste("qsub -l nodes=2,ncpus=10,mem=10gb,walltime=0:15:00 -q <queque.name> -o", 
    qsub.log.file,
    "-e", qsub.err.file , 
    cmd1)

        cat(qsub.cmd, "\n")
##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        
        system(qsub.cmd)
      }
      
    } else {
      ## the results folder for the current batch exists
      files <- list.files(path = path1)
      dt1   <- subset(dt_cnvr_raw, chr == chr1 & batch == batch1)
      dt1.failed <- subset(dt1, !fname %in% files)
      
      if (nrow(dt1.failed) == 0) {
        cat("CHR:", chr1, "BATCH:", batch1, "TOTAL:", freq1, "SUCCEED!\n")
      
      } else {  
        cat("CHR:", chr1, "BATCH:", batch1, "TOTAL:", freq1, "FAILED:", nrow(dt1.failed), "\n")
        
        if ( !dir.exists(paths = path_main_failed) ) {
          dir.create(path = path_main_failed, showWarnings = F, recursive = T)
        }
        
        fname.failed <- paste0("cnvrs_error_chr_", chr1, "_batch_", batch1, ".txt")
        write.table(data.frame(CNVR_ID = dt1.failed$CNVR_ID, stringsAsFactors = F),
                    file = file.path(path_main_failed, paste0("cnvrs_error_chr_", chr1, "_batch_", batch1, ".txt")),
                    col.names = T, row.names = F, quote = F)
        
        if (flag == 1) {
          cmd1 = paste(cmd, "--chr", chr1, "--batch", batch1, "--type", 1)
          
##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
## configure based on your system
    qsub.cmd <- paste("qsub -l nodes=2,ncpus=10,mem=10gb,walltime=0:15:00 -q <queque.name> -o", 
    qsub.log.file,
    "-e", qsub.err.file , 
    cmd1)
##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          
          cat(qsub.cmd, "\n")
          system(qsub.cmd)
        }
      }
    }
  }
}

# main runing function --------------------------------------------
check_jobs(path_main = path_main_pred, 
           dt_cnvr_raw = dt_cnvr_raw, 
           flag = flag, 
           path_main_failed = path_main_failed,
           path_joblog = path_joblog)



