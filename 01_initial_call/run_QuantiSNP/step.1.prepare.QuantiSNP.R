#!/hpc/packages/minerva-common/R/3.3.1/lib64/R/bin/Rscript --vanilla

suppressPackageStartupMessages(require(optparse))

option_list <- list(
  make_option(c("-i", "--input"), action = "store", default = NA, type = "character",
              help = "data folder for runing QuantiSNP"),
  make_option(c("-o", "--output"), action = "store", default = NA, type = "character",
              help = "output folder for QuantiSNP result")
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.na(opt$input) | is.na(opt$output)) {
  stop("All input and output arguments must be supplied.")
}

path_output <- opt$output
path_dat <- opt$input

## gender file
path_input <- "/sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/run.QuantiSNP/dat"
dat_gender <- readRDS(file = file.path(path_input, "Sample_ID_transform_detail.rds"))
dat_gender <- dat_gender[, c("Sample_ID_new", "Gender")]  ## use new sampleID gender and batch information
names(dat_gender) <- c("Sample_ID", "Gender")

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
  
  bsub.cmd <- paste("bsub -n 2 -W 02:00 -R 'rusage[mem=5000]' -P acc_haok01a", ##-R 'span[ptile=6]' 
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

