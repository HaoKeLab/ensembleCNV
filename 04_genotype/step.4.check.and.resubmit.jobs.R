#!/hpc/packages/minerva-common/R/3.3.1/lib64/R/bin/Rscript --vanilla

## check submited jobs until all jobs finished

suppressMessages(require(optparse))

option_list = list(
  make_option(c("-t", "--type"), action = "store", type = "integer", default = NA,
              help = "type 0 for only print CNVR runing state / type 1 for doing resubmit jobs for unfinished CNVR.\n")
)

opt = parse_args(OptionParser(option_list = option_list))
if ( is.na(opt$type) ) {
  stop("type only parameter must be supplied.( --help for detail)")
}

flag = as.integer(opt$type)  ## for 0 or 1

# resubmit jobs but exclude already finished CNVR
file_cnvr <- "cnvrs_annotated_batch.rds"  # with batch annotated
path_cnvr <- "/sc/orga/projects/haok01a/chengh04/Food_Allergy/FA_test/FA_1000/res/"


dt_cnvr_raw <- readRDS(file = file.path(path_cnvr, file_cnvr)) #
dt_cnvr_raw <- dt_cnvr_raw[order(dt_cnvr_raw$chr, dt_cnvr_raw$batch), ]
# add fname columns
dt_cnvr_raw$fname <- paste0(dt_cnvr_raw$CNVR_ID, "_pred.rds") 

tbl_raw <- table(dt_cnvr_raw$chr, dt_cnvr_raw$batch)
dt_freq_raw <- as.data.frame(tbl_raw)
names(dt_freq_raw) <- c("chr", "batch", "Freq")

dt_freq_raw <- subset(dt_freq_raw, Freq != 0)  ## subset no 0 freq part
dt_freq_raw <- dt_freq_raw[order(dt_freq_raw$chr, dt_freq_raw$batch), ]

## main predict path
path_main_pred <- "/sc/orga/projects/haok01a/chengh04/Food_Allergy/FA_test/FA_1000/res/pred"
path_main_failed <- "/sc/orga/projects/haok01a/chengh04/Food_Allergy/FA_test/FA_1000/res/cnvrs_failed"

# function ----------------------------------------------------------
# submit jobs
submit_jobs <- function(chr1, batch1, type1) {
  
  path_code <- "/sc/orga/projects/haok01a/chengh04/Food_Allergy/FA_test/code/regenotype"
  name_script <- "step.2.runme.minerva.each.chr.each.batch.R"
  
  script <- file.path(path_code, name_script)
  
  cmd1 = paste(script, "-c", chr1, "-b", batch1, "-t", 0)
  bsub.cmd <- paste("bsub -n 2 -W 05:00 -R 'rusage[mem=10000]' -P acc_haok01a",
                    "-q premium", shQuote(cmd1))
  
  cat(cmd1, "\n")  ## print out
  system(bsub.cmd)
  Sys.sleep(0.1)
}

# check if all CNVR finished --------------------------------------------------------------
check_jobs <- function(path_main, dt_cnvr_raw, flag, path_main_failed) {
  
  ## remove folder and create
  system(paste("rm -rf", path_main_failed)) ## remove all previous results
  
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
        submit_jobs(chr1 = chr1, batch1 = batch1, type1 = 0) # do not exist folder
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
          dir.create(path = path_main_failed, showWarnings = FALSE)
        }
        
        fname.failed <- paste0("cnvrs_chr_", chr1, "_batch_", batch1, "_failed.rds")
        saveRDS(dt1.failed, file = file.path(path_main_failed, fname.failed))
        
        if (flag == 1) {
          submit_jobs(chr1 = chr1, batch1 = batch1, type1 = 1)
        }
        
      }
      
    }
    
  }
  
}

# main run function --------------------------------------------
check_jobs(path_main = path_main_pred, dt_cnvr_raw = dt_cnvr_raw, 
           flag = flag, path_main_failed = path_main_failed)



