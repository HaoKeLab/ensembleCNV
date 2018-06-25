#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(optparse))

## function ------------------------------------------------------------------
run.quantisnp <- function(path_output, path_dat, sample_name, gender) {
  
  ## define program variables
  EMITERS    <- "10"        ## number of EM iterations to use during training
  LSETTING   <- "2000000"   ## characteristic CNV length parameter
  GCDIR      <- "/hpc/packages/minerva-common/quantisnp/data/b37/"   ## set path to GC data files (contents of gc_data.zip)
  PARAMSFILE <- "/hpc/packages/minerva-common/quantisnp/2.3/quantisnp/config/params.dat"      ## path to parameters file
  LEVELSFILE <- "/hpc/packages/minerva-common/quantisnp/2.3/quantisnp/config/levels-hd.dat"   ## path to levels file
  CHRRANGE   <- "1:23"   ## path to parameters file
  CHRX       <- "23"     ## which chromosome is X?
  OUTDIR     <- file.path(path_output, sample_name)    ## output directory
  SAMPLEID   <- sample_name ## sample name
  GENDER     <- gender      ## sample gender
  INFILE     <- file.path(path_dat, paste0(sample_name, ".txt"))   ## input data file
  MCRROOT    <- "/hpc/packages/minerva-common/quantisnp/MATLAB_RT/lib/v79/"   ## set path to MCR Run-Time Libraries
  
  if(!file.exists(OUTDIR)) dir.create(OUTDIR)
  
  cmd <- paste("/hpc/packages/minerva-common/quantisnp/2.3/quantisnp/linux64/run_quantisnp2.sh",
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
  
  bsub.cmd <- paste("bsub -n 2 -W 10:00 -R 'rusage[mem=5000]' -P acc_haok01a", ##-R 'span[ptile=6]' 
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
              help = "data path."),
  make_option(c("-r", "--result"), default = NA, type = "character", action = "store",
              help = "call CNV results path.")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.na(opt$data) | is.na(opt$result)) {
  stop("Three input argument must be supplied.")
}

# get paras
path_data <- opt$data
path_res <- opt$result


## check QuantiSNP jobs failed or not finished
# path_input <- "/sc/orga/projects/haok01a/chengh04/FMD.GWAS/genotype_QC/2017-recluster-QC/REPORT" # for 2017
# path_input <- "/sc/orga/projects/haok01a/chengh04/FMD.GWAS/genotype_QC/2015-recluster-QC/REPORT" ## for 2015
# name_input <- "QC.samples.xls"

# dat_gender <- read.table(file = file.path(path_input, name_input),
#                          sep = "\t", header = TRUE, check.names = FALSE, as.is = TRUE)

# dat_gender <- read.table(file = file_gender,
#                          sep = "\t", header = TRUE, check.names = FALSE, as.is = TRUE)
# 
# 
# samples <- dat_gender$IID  ## all folder names
# genders <- dat_gender$sex.inferred ## all genders

## gender file
path_input <- "/sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/run.QuantiSNP/dat"
dat_gender <- readRDS(file = file.path(path_input, "Sample_ID_transform_detail.rds"))
dat_gender <- dat_gender[, c("Sample_ID_new", "Gender")]
names(dat_gender) <- c("Sample_ID", "Gender")

samples <- dat_gender$Sample_ID
genders <- dat_gender$Gender

# for 2017
# path_res <- "/sc/orga/projects/haok01a/chengh04/FMD.GWAS/callCNV/Kovacic_128samples_042717/QuantiSNP/res"
# path_dat <- "/sc/orga/projects/haok01a/chengh04/FMD.GWAS/callCNV/Kovacic_128samples_042717/QuantiSNP/data"
# for 2015
# path_res <- "/sc/orga/projects/haok01a/chengh04/FMD.GWAS/callCNV/Valentina_112samples_120214/QuantiSNP/res"
# path_dat <- "/sc/orga/projects/haok01a/chengh04/FMD.GWAS/callCNV/Valentina_112samples_120214/QuantiSNP/data"

for (i in 1:length(samples)) {
  
  sample_name <- samples[i]
  gender <- genders[i]
  path_sample1 <- file.path(path_res, sample_name)
  
  if (dir.exists(paths = path_sample1)) {
    
    # check if have generated .cnv file
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