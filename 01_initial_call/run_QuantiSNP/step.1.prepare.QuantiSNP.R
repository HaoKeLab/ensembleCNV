#!/usr/bin/env Rscript

## The script was used to run QuantiSNP on Minerva high performance cluster.
## You need to modifiy it according to the system you are using if you would like to use it.
## Please refer to original QuantiSNP documents (https://sites.google.com/site/quantisnp/) for more information 

path_to_quantisnp <- ""  ## where QuantiSNP is installed

suppressPackageStartupMessages(require(optparse))

option_list <- list(
  make_option(c("-i", "--input"), action = "store", default = NA, type = "character",
              help = "data folder for runing QuantiSNP"),
  make_option(c("-o", "--output"), action = "store", default = NA, type = "character",
              help = "output folder for QuantiSNP results")
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.na(opt$input) | is.na(opt$output)) {
  stop("All input and output arguments must be supplied.")
}

path_dat <- opt$input
path_output <- opt$output

## path to other auxiliary information (e.g. gender file)
path_input <- "" 

## gender file: in tab-delimited format and has two columns: Sample_ID and Gender
## for example
# Sample_ID	Gender
# sample_1	female
# sample_2	male

dat_gender <- read.delim(file = file.path(path_input, "gender_file.txt"), as.is = TRUE)

cat("rows of dat_gender:", nrow(dat_gender), "\n") ## number of samples

for (i in 1:nrow(dat_gender)) {
  
  sample_name <- as.character(dat_gender$Sample_ID[i])
  gender <- tolower(as.character(dat_gender$Gender[i]))
  
  # check 
  res_files <- list.files(path = file.path(path_output, sample_name))
  idx <- grep(pattern = "cnv", res_files)
  if (length(idx) >0) {
    cat("i:", sample_name, "\n")
    next
  }
  
  ## define program variables
  EMITERS    <- "10"        ## number of EM iterations to use during training
  LSETTING   <- "2000000"   ## characteristic CNV length parameter
  GCDIR      <- file.path(path_to_quantisnp, "data/b37/")              ## path to GC data files (contents of gc_data.zip)
  PARAMSFILE <- file.path(path_to_quantisnp, "config/params.dat")      ## path to parameters file
  LEVELSFILE <- file.path(path_to_quantisnp, "config/levels-hd.dat")   ## path to levels file
  MCRROOT    <- file.path(path_to_quantisnp, "MATLAB_RT/lib/v79/")     ## path to MCR Run-Time Libraries
  CHRRANGE   <- "1:23"   ## chromosomes
  CHRX       <- "23"     ## which chromosome is X?
  OUTDIR     <- file.path(path_output, sample_name)    ## output directory
  SAMPLEID   <- sample_name ## sample name
  GENDER     <- gender      ## sample gender
  INFILE     <- file.path(path_dat, paste0(sample_name, ".txt"))   ## input data file generated with finalreport_to_QuantiSNP.pl

  
  if (!file.exists(OUTDIR)) dir.create(OUTDIR)
  
  cmd <- paste(file.path(path_to_quantisnp, "linux64/run_quantisnp.sh"),
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
  
  bsub.cmd <- paste("bsub -n 2 -W 02:00 -R 'rusage[mem=5000]' -P [account]",
                    "-J", job.name,
                    "-q premium",
                    "-oo", log.file,
                    "-eo", err.file ,
                    shQuote(cmd))
  
  cat("i =", i, bsub.cmd, "\n")
  system(bsub.cmd)
  Sys.sleep(0.1)
  
  cat("i = ", i , sample_name, "\n")
}

