#!/usr/bin/env Rscript

## The script was used to run QuantiSNP on Minerva high performance cluster.
## You need to modifiy it according to the system you are using if you would like to use it.
## Please refer to original QuantiSNP documents (https://sites.google.com/site/quantisnp/) for more information 

suppressPackageStartupMessages(require(optparse))

## function ------------------------------------------------------------------
run.quantisnp <- function(path_output, path_dat, sample_name, gender) {
  
  ## define program variables
  EMITERS    <- "10"        ## number of EM iterations to use during training
  LSETTING   <- "2000000"   ## characteristic CNV length parameter
  GCDIR      <- "path_to_quantisnp/data/b37/"   ## set path to GC data files (contents of gc_data.zip)
  PARAMSFILE <- "path_to_quantisnp/config/params.dat"      ## path to parameters file
  LEVELSFILE <- "path_to_quantisnp/config/levels-hd.dat"   ## path to levels file
  CHRRANGE   <- "1:23"   ## chromosome
  CHRX       <- "23"     ## which chromosome is X?
  OUTDIR     <- file.path(path_output, sample_name)    ## output directory
  SAMPLEID   <- sample_name ## sample name
  GENDER     <- gender      ## sample gender
  INFILE     <- file.path(path_dat, paste0(sample_name, ".txt"))   ## input data file
  MCRROOT    <- "path_to_quantisnp/MATLAB_RT/lib/v79/"   ## set path to MCR Run-Time Libraries
  
  if (!file.exists(OUTDIR)) dir.create(OUTDIR)
  
  cmd <- paste("path_to_quantisnp/linux64/run_quantisnp.sh",
               MCRROOT, 
               paste("--chr", CHRRANGE),
               paste("--outdir", OUTDIR), 
               paste("--sampleid", SAMPLEID),
               paste("--gender", GENDER), 
               paste("--emiters", EMITERS), 
               paste("--lsetting", LSETTING), 
               paste("--gcdir", GCDIR),
               "--plot", 
               "--genotype", 
               paste("--config", PARAMSFILE), 
               paste("--levels", LEVELSFILE), 
               paste("--input-files", INFILE), 
               paste("--chrX", CHRX), 
               "--doXcorrect")
  
  job.name <- sample_name
  log.file <- file.path(OUTDIR, paste0(sample_name, ".quantisnp.log"))
  err.file <- file.path(OUTDIR, paste0(sample_name, ".quantisnp.err"))
  
  bsub.cmd <- paste("bsub -n 2 -W 10:00 -R 'rusage[mem=5000]' -P [account]", 
                    "-J", job.name,
                    "-q premium",
                    "-oo", log.file,
                    "-eo", err.file ,
                    shQuote(cmd))
  
  cat(bsub.cmd, "\n")
  system(bsub.cmd)
  
}

## ===============================================================================================

option_list <- list(
  make_option(c("-d", "--data"), default = NA, type = "character", action = "store",
              help = "data folder for runing QuantiSNP."),
  make_option(c("-r", "--result"), default = NA, type = "character", action = "store",
              help = "path to CNV results generated in the first step.")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.na(opt$data) | is.na(opt$result)) {
  stop("Three input argument must be supplied.")
}

# get paras
path_data <- opt$data
path_res <- opt$result

## gender file
# change path_input and filename here 
# file in .rds format, must have two columns: Sample_ID and Gender
path_input <- "" ## path to other auxiliary information
dat_gender <- readRDS(file = file.path(path_input, "gender_file.rds"))

cat("rows of dat_gender:", nrow(dat_gender), "\n") ## number of samples

samples <- dat_gender$Sample_ID
genders <- dat_gender$Gender

for (i in 1:length(samples)) {
  
  sample_name <- samples[i]
  gender <- genders[i]
  path_sample1 <- file.path(path_res, sample_name)
  
  if (dir.exists(paths = path_sample1)) {
    
    # check if .cnv file have been generated
    files <- list.files(path = path_sample1)
    idx1 <- grep(pattern = ".cnv", files)
    if (length(idx1) == 1) {
      cat("Sample_ID:", sample_name, "SUCCESS.\n")
    } else {
      cat("Sample_ID:", sample_name, "FAILED.\n")
      run.quantisnp(path_output = path_res, path_dat = path_data, sample_name = sample_name, gender = gender)
    }

  } else {
    cat("Sample_ID:", sample_name, "FAILED.\n")
    run.quantisnp(path_output = path_res, path_dat = path_data, sample_name = sample_name, gender = gender)
  }
  
}