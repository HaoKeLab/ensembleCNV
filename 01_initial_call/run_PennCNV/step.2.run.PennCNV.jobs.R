#!/usr/bin/env Rscirpt

## The script was used to run PennCNV on Minerva high performance cluster.
## You need to modifiy it according to the system you are using if you would like to use it.
## Please refer to original PennCNV documents (http://penncnv.openbioinformatics.org/en/latest/) for more information 

suppressMessages({
  require( optparse, quietly = TRUE)
})

option_list <- list(
  make_option(c("-a", "--data"), action = "store", default = NA, type = "character",
              help = "path to text data file for each sample."),
  make_option(c("-d", "--main"), action = "store", default = NA, type = "character",
              help = "main path for each sample's list file and res."),
  make_option(c("-f", "--pfb"), action = "store", default = NA, type = "character",
              help = "pfb file."),
  make_option(c("-g", "--gcmodel"), action = "store", default = NA, type = "character",
              help = "gcmodel file."),
  make_option(c("-m", "--hmm"), action = "store", default = NA, type = "character",
              help = "HMM model file.")
)

opt = parse_args(OptionParser(option_list = option_list))

path_dat     <- opt$dat
path_main    <- opt$main
file_pfb     <- opt$pfb
file_gcmodel <- opt$gcmodel
file_hmm     <- opt$hmm

if (any(is.na(c(path_dat, path_main, file_pfb, file_gcmodel, file_hmm)))) {
  stop("all parameters must be supplied. (--help for details)")
}

# create path -------------------------------------------------------------

path_list  <- file.path(path_main, "list")
path_res   <- file.path(path_main, "res")

if ( !dir.exists(path_list) ) {
  dir.create(path = path_list, showWarnings = FALSE, recursive = TRUE)
}
if ( !dir.exists(path_res) ) {
  dir.create(path = path_res, showWarnings = FALSE, recursive = TRUE)
}


# generate list.txt for each sample ---------------------------------------

sample_files <- list.files(path = path_dat)

cat("number of samples:", length(sample_files), "\n")

for ( i in 1:length(sample_files) ) {

  sample_file <- sample_files[i]
  dat1 <- data.frame(file_name = file.path(path_dat, sample_file), ## add whole path information
                     stringsAsFactors = FALSE)
  write.table(dat1, file = file.path(path_list, sample_file), sep = "\t",
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}


# cmd_PennCNV -------------------------------------------------------------

cmd_PennCNV <- function(file_hmm, file_pfb, file_gcmodel,
                        filename_sample, path_list, path_res_sample) {

  file_list <- file.path(path_list, filename_sample)

  samplename <- gsub(pattern = ".txt$", replacement = "", filename_sample)
  file_log   <- file.path(path_res_sample, paste0(samplename, ".log"))
  file_rawcnv <- file.path(path_res_sample, paste0(samplename, ".rawcnv"))

  cmd <- paste("/hpc/packages/minerva-common/penncnv/2011Jun16/bin/detect_cnv.pl",
               "-test --confidence",
               "-hmm", file_hmm,
               "-pfb", file_pfb,
               "-gcmodel", file_gcmodel,
               "-list", file_list,
               "-log", file_log,
               "-out", file_rawcnv)
  cmd
}

cmd_submitjob <- function(cmd.sample, samplename) {

  bsub.cmd <- paste("bsub -n 2 -W 00:30 -R 'rusage[mem=5000]' -P acc_haok01a", ##-R 'span[ptile=6]'
                    "-J", samplename,
                    "-q premium",
                    shQuote(cmd.sample))
  bsub.cmd
}

# main loop ---------------------------------------------------------------

for ( i in 1:length(sample_files) ) {

  sample_file <- sample_files[i]
  samplename <- gsub(pattern = ".txt$", replacement = "", sample_file)

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

  system(cmd.job)
  Sys.sleep(0.1)

}


