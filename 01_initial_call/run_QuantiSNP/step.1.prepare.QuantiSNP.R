#!/usr/bin/env Rscript

## NOTE: The scripts embraced by "##<<<... ##>>>..." need to be specified based on your system

## The script was used to run QuantiSNP on Minerva high performance cluster.
## You need to modifiy it according to the system you are using if you would like to use it.
## Please refer to original QuantiSNP documents (https://sites.google.com/site/quantisnp/) for more information 

## sample file: in tab-delimited format and has two columns: Sample_ID and Gender
## for example
# Sample_ID	Gender
# sample_1	Female
# sample_2	Male

suppressPackageStartupMessages(require(optparse))

option_list <- list(
  make_option(c("-q", "--quantisnp"), action = "store", default = NA, type = "character",
              help = "path to QuantiSNP installation folder."),  
  make_option(c("-i", "--input"), action = "store", default = NA, type = "character",
              help = "data folder for runing QuantiSNP"),
  make_option(c("-s", "--sample"), action = "store", default = NA, type = "character",
              help = "sample file with Sample_ID and Gender information for runing QuantiSNP"),            
  make_option(c("-o", "--output"), action = "store", default = NA, type = "character",
              help = "output folder for QuantiSNP results")
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.na(opt$input) | is.na(opt$output)) {
  stop("All input and output arguments must be supplied.")
}

path_quantisnp <- opt$quantisnp
path_dat       <- opt$input
sample_file    <- opt$sample
path_output    <- opt$output

dat_sample <- read.delim(file = sample_file, as.is = TRUE)

cat("number of rows of sample table:", nrow(dat_sample), "\n") ## number of samples

for (i in 1:nrow(dat_sample)) {
  
  sample_name <- as.character(dat_sample$Sample_ID[i])
  gender <- tolower(as.character(dat_sample$Gender[i]))
  ## must change Female => female and Male => male
  
  ## check if folder exists
  res_files <- list.files(path = file.path(path_output, sample_name))
  idx <- grep(pattern = "cnv", res_files)
  if (length(idx) >0) {
    cat("i:", sample_name, "\n")
    next
  }
  
  ## define program variables
  EMITERS    <- "10"        ## number of EM iterations to use during training
  LSETTING   <- "2000000"   ## characteristic CNV length parameter
##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   
  GCDIR      <- file.path(path_quantisnp, "data/b37/")                        ## path to GC data files (contents of gc_data.zip)
  PARAMSFILE <- file.path(path_quantisnp, "quantisnp/config/params.dat")      ## path to parameters file
  LEVELSFILE <- file.path(path_quantisnp, "quantisnp/config/levels-hd.dat")   ## path to levels file
  MCRROOT    <- file.path(path_quantisnp, "v79/")                             ## path to MCR Run-Time Libraries
##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  
  CHRRANGE   <- "1:23"   ## chromosomes
  CHRX       <- "23"     ## which chromosome is X?
  OUTDIR     <- file.path(path_output, sample_name)    ## output directory
  SAMPLEID   <- sample_name ## sample name
  GENDER     <- gender      ## sample gender
  INFILE     <- file.path(path_dat, paste0(sample_name, ".txt"))   ## input data file generated with finalreport_to_QuantiSNP.pl

  
  if (!file.exists(OUTDIR)) dir.create(OUTDIR)

##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  
  cmd <- paste(file.path(path_quantisnp, "quantisnp/linux64/run_quantisnp2.sh"),
##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
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

##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
## configure based on your system  
  bsub.cmd <- paste("bsub -n 2 -W 02:00 -R 'rusage[mem=5000]' -P <account>",
                    "-J", job.name,
                    "-q premium",
                    "-oo", log.file,
                    "-eo", err.file ,
                    shQuote(cmd))
##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  cat("i =", i, bsub.cmd, "\n")
  system(bsub.cmd)
  Sys.sleep(0.1)
  
  cat("i = ", i , sample_name, "\n")
}

