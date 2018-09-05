#!/usr/bin/env Rscript

suppressMessages(require(optparse))

option_list = list(
  make_option(c("-p", "--datapath"), action = "store", type = "character", default = NA,
              help = "path contain all running needed data"),
  make_option(c("-o", "--resultpath"), action = "store", type = "character", default = NA,
              help = "path save all results"),
  make_option(c("-m", "--matrixpath"), action = "store", type = "character", default = NA,
              help = "path save all LRR and BAF chromosome based data"),
  make_option(c("-s", "--sourcefile"), action = "store", type = "character", default = NA,
              help = "path contain all scripts need to be soucred"),
  make_option(c("-d", "--duplicates"), action = "store_true", default = FALSE,
              help = "file of duplicate pairs."),
  make_option(c("-n", "--plot"), action = "store_true", default = FALSE,
              help = "plot png or not"),
  make_option(c("-r", "--script"), action = "store", type = "character", default = NA,
              help = "script path named as CNV.genotype.one.chr.one.batch.HC.R."),
  make_option(c("-l", "--joblog"), action = "store", type = "character", default = NA,
              help = "path save all jobs log."),
  make_option(c("-f", "--flag"), action = "store", type = "integer", default = NA,
              help = "0 for only print CNVR runing state/1 for doing resubmit jobs for unfinished CNVR ")
)



opt = parse_args(OptionParser(option_list = option_list))
pars = c(opt$type, opt$datapath, opt$resultpath, opt$joblog,
         opt$matrixpath, opt$sourcefile, opt$script, opt$flag)

if ( any(is.na(pars)) ) {
  stop("All parameters must be supplied. (--help for detail)")
}

flag = as.integer( opt$flag )  ## for 0 or 1

# resubmit jobs but exclude already finished CNVR
file_cnvr <- "cnvr_batch.txt"  ## with batch information
dt_cnvr_raw <- read.delim(file = file.path(opt$datapath, file_cnvr), as.is = TRUE)
# dt_cnvr_raw <- readRDS(file = file.path(path_cnvr, file_cnvr)) #
dt_cnvr_raw <- dt_cnvr_raw[order(dt_cnvr_raw$chr, dt_cnvr_raw$batch), ]
# add fname columns
dt_cnvr_raw$fname <- paste0(dt_cnvr_raw$CNVR_ID, "_pred.rds") 

tbl_raw <- table(dt_cnvr_raw$chr, dt_cnvr_raw$batch)
dt_freq_raw <- as.data.frame(tbl_raw)
names(dt_freq_raw) <- c("chr", "batch", "Freq")

dt_freq_raw <- subset(dt_freq_raw, Freq != 0)  ## subset no 0 freq part
dt_freq_raw <- dt_freq_raw[order(dt_freq_raw$chr, dt_freq_raw$batch), ]

## main predict path
# path_main_pred <- "/sc/orga/projects/haok01a/chengh04/Food_Allergy/FA_test/FA_1000/res/pred"
# path_main_failed <- "/sc/orga/projects/haok01a/chengh04/Food_Allergy/FA_test/FA_1000/res/cnvrs_failed"

path_main_pred <- file.path(opt$resultpath, "pred")
path_main_failed <- file.path(opt$resultpath, "cnvrs_error")

# create script
script <- file.path(opt$script, "CNV.genotype.one.chr.one.batch.HC.R")
cmd    <- paste("Rscript", script, 
                "--datapath", opt$datapath,
                "--resultpath", opt$resultpath,
                "--matrixpath", opt$matrixpath,
                "--sourcefile", opt$sourcefile)

if ( opt$duplicates ) cmd <- paste(cmd, "--duplicates")
if ( opt$png ) cmd <- paste(cmd, "--png")

# function ----------------------------------------------------------
# submit jobs
# submit_jobs <- function(chr1, batch1, type1) {
#   
#   path_code <- "/sc/orga/projects/haok01a/chengh04/Food_Allergy/FA_test/code/regenotype"
#   name_script <- "step.2.runme.minerva.each.chr.each.batch.R"
#   
#   script <- file.path(path_code, name_script)
#   
#   cmd1 = paste(script, "-c", chr1, "-b", batch1, "-t", 0)
#   bsub.cmd <- paste("bsub -n 2 -W 05:00 -R 'rusage[mem=10000]' -P acc_haok01a",
#                     "-q premium", shQuote(cmd1))
#   
#   cat(cmd1, "\n")  ## print out
#   system(bsub.cmd)
#   Sys.sleep(0.1)
# }

# check if all CNVR finished --------------------------------------------------------------
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
    freq1 <- dt_freq_raw$Freq[i]  # total freq 
    
    foldername1 <- paste0("chr_", chr1, "_batch_", batch1)
    path1 <- file.path(path_main, foldername1)  # folder name
    
    if ( !dir.exists(paths = path1) ) { 
      cat("CHR:", chr1, "BATCH:", batch1, "Failed jobs must be resubmit.\n")
      # submit jobs
      if (flag == 1) {
        
        cmd1 = paste(cmd, "--chr", chr1, "--batch", batch1, "--type", 0)
        bsub.cmd = paste("bsub -n 2 -W 10:00 -R 'rusage[mem=20000]' -P [account]", ## need to modify based on specific system
                         "-e", file.path(path_job_error, paste0("chr_", chr1, "_batch_", batch1, ".e")), 
                         "-o", file.path(path_job_out, paste0("chr_", chr1, "_batch_", batch1, ".o")),
                         "-q premium", shQuote(cmd1))
        cat(bsub.cmd, "\n")
        system(bsub.cmd)
        
      }
      
    } else {
      
      files <- list.files(path = path1)
      dt1   <- subset(dt_cnvr_raw, chr == chr1 & batch == batch1)
      dt1.failed <- subset(dt1, !fname %in% files)
      
      if (nrow(dt1.failed) == 0) {
        cat("CHR:", chr1, "BATCH:", batch1, "TOTAL:", freq1, "SUCCESSED!!!!!!!\n")
      } else {
        
        cat("CHR:", chr1, "BATCH:", batch1, "TOTAL:", freq1, "FAILED:", nrow(dt1.failed), "\n")
        
        if ( !dir.exists(paths = path_main_failed) ) {
          dir.create(path = path_main_failed, showWarnings = F, recursive = T)
        }
        
        # fname.failed <- paste0("cnvrs_chr_", chr1, "_batch_", batch1, "_failed.rds")
        # saveRDS(dt1.failed, file = file.path(path_main_failed, fname.failed))
        fname.failed <- paste0("cnvrs_error_chr_", chr1, "_batch_", batch1, ".txt")
        write.table(data.frame(CNVR_ID = dt1.failed$CNVR_ID, stringsAsFactors = F),
                    file = file.path(path_main_failed, paste0("cnvrs_error_chr_", chr1, "_batch_", batch1, ".txt")),
                    col.names = T, row.names = F, quote = F)
        
        if (flag == 1) {
          cmd1 = paste(cmd, "--chr", chr1, "--batch", batch1, "--type", 1)
          bsub.cmd = paste("bsub -n 2 -W 10:00 -R 'rusage[mem=20000]' -P [account]", ## need to modify based on specific system
                           "-e", file.path(path_job_error, paste0("chr_", chr1, "_batch_", batch1, ".e")), 
                           "-o", file.path(path_job_out, paste0("chr_", chr1, "_batch_", batch1, ".o")),
                           "-q premium", shQuote(cmd1))
          cat(bsub.cmd, "\n")
          system(bsub.cmd)
        }
      }
    }
  }
}

# main runing function --------------------------------------------
check_jobs(path_main = path_main_pred, dt_cnvr_raw = dt_cnvr_raw, 
           flag = flag, path_main_failed = path_main_failed,
           path_joblog = path_joblog)



