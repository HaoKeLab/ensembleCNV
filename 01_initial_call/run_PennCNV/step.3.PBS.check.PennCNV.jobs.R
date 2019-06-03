#!/usr/bin/env Rscript

##############################
###### THIS IS MODIFIED ###### 
##############################

## NOTE: here I reorganized the sript in order to makes it work with PBS and the qsub command --S.M.
## All the change should be commented

## NOTE: The scripts embraced by "##<<<... ##>>>..." need to be specified based on your system

## The script was used to run PennCNV on Minerva high performance cluster.
## You need to modifiy it according to the system you are using if you would like to use it.
## Please refer to original PennCNV documents (http://penncnv.openbioinformatics.org/en/latest/) for more information

suppressMessages({
  require( optparse, quietly = TRUE)
})

options(warn = 2)

option_list <- list(
  make_option(c("-p", "--penncnv"), action = "store", default = NA, type = "character",
              help = "path to PennCNV installation folder."),  
  make_option(c("-a", "--data"), action = "store", default = NA, type = "character",
              help = "path to tab-delimit text data files for each sample."),
  make_option(c("-d", "--wkdir"), action = "store", default = NA, type = "character",
              help = "working directory."),
  make_option(c("-f", "--pfb"), action = "store", default = NA, type = "character",
              help = "pfb file."),
  make_option(c("-g", "--gcmodel"), action = "store", default = NA, type = "character",
              help = "gcmodel file."),
  make_option(c("-m", "--hmm"), action = "store", default = NA, type = "character",
              help = "HMM model file.")
)

opt = parse_args(OptionParser(option_list = option_list))

path_penncnv <- opt$penncnv
path_data    <- opt$data
path_wkdir   <- opt$wkdir
file_pfb     <- opt$pfb
file_gcmodel <- opt$gcmodel
file_hmm     <- opt$hmm

if (any(is.na(c(path_data, path_wkdir, file_pfb, file_gcmodel, file_hmm)))) {
  stop("All parameters must be supplied. (--help for details)")
}

path_list  <- file.path(path_wkdir, "list")
path_res  <- file.path(path_wkdir, "res")  ## PennCNV results folder

# submit jobs functions ---------------------------------------------------

cmd_PennCNV <- function(file_hmm, file_pfb, file_gcmodel, 
                        filename_sample, path_list, path_res_sample) {
  
  file_list <- file.path(path_list, sub("\\.txt$", ".list", filename_sample))
  
  samplename <- gsub(pattern = "\\.txt$", replacement = "", filename_sample)
  file_log   <- file.path(path_res_sample, paste0(samplename, ".log"))
  file_rawcnv <- file.path(path_res_sample, paste0(samplename, ".rawcnv"))
  
   cmd <- paste("-v HMM=", file_hmm,
                ",PFB=", file_pfb,
                ",GC=", file_gcmodel,
                ",LIST=", file_list,
                ",LOG=", file_log,
                ",OUT=", file_rawcnv,
		",PENNCNV=", path_penncnv,
	        " ",
                file.path(path_penncnv, "bin/run_detect_cnv.sh"), 
	        sep="")
  
  cmd
} 

cmd_submitjob <- function(cmd.sample, samplename) {

##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
## configure based on your system

	qsub.cmd <- paste("qsub -l nodes=2,ncpus=10,mem=10gb,walltime=0:15:00 -q <queque.name> -N", samplename, cmd.sample)

##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  qsub.cmd
}

# main loop ---------------------------------------------------------------

sample_files <- list.files(path = path_data)
cat("number of samples:", length(sample_files), "\n")

n.success <- 0
n.fail <- 0
for ( i in 1:length(sample_files) ) {
  
  sample_file <- sample_files[i]
  samplename  <- gsub(pattern = "\\.txt$", replacement = "", sample_file)
  
  path_res_sample <- file.path(path_res, samplename)
  file_rawcnv <- file.path(path_res_sample, paste0(samplename, ".rawcnv"))
  
  flag.folder <- dir.exists(paths = path_res_sample)
  flag.rawcnv <- file.exists(file_rawcnv)
  
  if ( flag.folder & flag.rawcnv ) {
    cat("Sample_ID:", samplename, "SUCCESS\n")
    n.success <- n.success + 1
  } else {
    
    cat("Sample_ID:", samplename, "FAILED\n")
    dir.create(path = path_res_sample, showWarnings = FALSE, recursive = TRUE)
    
    cmd.sample <- cmd_PennCNV(file_hmm = file_hmm,
                              file_pfb = file_pfb,
                              file_gcmodel = file_gcmodel,
                              filename_sample = sample_file,
                              path_list = path_list,
                              path_res_sample = path_res_sample)
    
    cmd.job   <- cmd_submitjob(cmd.sample = cmd.sample, samplename = samplename)
   cat(cmd.job, "\n")
 
    system(cmd.job)
    Sys.sleep(0.1)
    
    n.fail <- n.fail + 1
    
  }
}

cat("total number of samples:", length(sample_files),
    "number of success:", n.success,
    "number of fail:", n.fail, "\n")


