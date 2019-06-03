#!/usr/bin/env Rscirpt

##############################
###### THIS IS MODIFIED ###### 
##############################

## NOTE: here I reorganized the sript in order to makes it work with PBS and the qsub command --S.M.
## All the change should be commented

## NOTE: that it is necessary to copy the script run_detect_cnv.sh into $PENNCNV/bin/ to make this work 

## NOTE: The scripts embraced by "##<<<... ##>>>..." need to be specified based on your system

## The script was used to run PennCNV on Minerva high performance cluster.
## You need to modifiy it according to the system you are using if you would like to use it.
## Please refer to original PennCNV documents (http://penncnv.openbioinformatics.org/en/latest/) for more information 


suppressMessages({
  require( optparse, quietly = TRUE)
})

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

## CHECK if the additional script has been placed in the correct location 
if (!file.exists(file.path(path_penncnv, "bin/run_detect_cnv.sh"))){
        stop("Place the additional script in the correct location ($PENNCNV/bin/run_detect_cnv.sh)")
        }


# create path -------------------------------------------------------------

path_list  <- file.path(path_wkdir, "list")
path_res   <- file.path(path_wkdir, "res") ## PennCNV raw results folder

if ( !dir.exists(path_list) ) {
  dir.create(path = path_list, showWarnings = FALSE, recursive = TRUE)
}
if ( !dir.exists(path_res) ) {
  dir.create(path = path_res, showWarnings = FALSE, recursive = TRUE)
}


# generate list.txt for each sample ---------------------------------------

sample_files <- list.files(path = path_data)

cat("number of samples:", length(sample_files), "\n")

for ( i in 1:length(sample_files) ) {

  sample_file <- sample_files[i]
  sample_list <- sub("\\.txt$", ".list", sample_file)
  
  dat1 <- data.frame(file_name = file.path(path_data, sample_file), ## add whole path information
                     stringsAsFactors = FALSE)
  write.table(dat1, file = file.path(path_list, sample_list), sep = "\t",
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}


# cmd_PennCNV -------------------------------------------------------------

cmd_PennCNV <- function(file_hmm, file_pfb, file_gcmodel,
                        filename_sample, path_list, path_res_sample) {

  file_list <- file.path(path_list, sub("\\.txt$", ".list", filename_sample))

  samplename <- gsub(pattern = "\\.txt$", replacement = "", filename_sample)
  
  file_log   <- file.path(path_res_sample, paste0(samplename, ".log"))
  file_rawcnv <- file.path(path_res_sample, paste0(samplename, ".rawcnv"))
	
## here main changes

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

## ADDED qsub.log and qsub.err files here
	qsub.log.file <- file.path(path_res_sample, paste0(samplename, ".log.qsub"))
        qsub.err.file <- file.path(path_res_sample, paste0(samplename, ".err.qsub"))

##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
## configure based on your system

qsub.cmd <- paste("qsub -l nodes=2,ncpus=10,mem=10gb,walltime=0:15:00 -q <queque.name> -o", qsub.log.file,"-e", qsub.err.file , "-N", samplename, cmd.sample)


##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  
  qsub.cmd
}

# main loop ---------------------------------------------------------------

for ( i in 1:length(sample_files) ) {

  sample_file <- sample_files[i]
  samplename <- gsub(pattern = "\\.txt$", replacement = "", sample_file)

  path_res_sample <- file.path(path_res, samplename)
  dir.create(path = path_res_sample, showWarnings = FALSE, recursive = TRUE)

  cat("Sample_ID:", samplename, "\n")

  cmd.sample <- cmd_PennCNV(file_hmm = file_hmm,
                            file_pfb = file_pfb,
                            file_gcmodel = file_gcmodel,
                            filename_sample = sample_file,
                            path_list = path_list,
                            path_res_sample = path_res_sample)


  cmd.job   <- cmd_submitjob(cmd.sample = cmd.sample, samplename = samplename)
  
  cat("i =", i, cmd.job, "\n")

  system(cmd.job)
  Sys.sleep(0.1)

}


