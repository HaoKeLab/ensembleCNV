#!/hpc/packages/minerva-common/R/3.3.1/lib64/R/bin/Rscript --vanilla

suppressPackageStartupMessages(require(optparse))

option_list <- list(
  make_option(c("-a", "--dat"), action = "store", default = NA, type = "character",
              help = "dat path from each sample for running PennCNV."),
  make_option(c("-b", "--main"), action = "store", default = NA, type = "character",
              help = "main path for each sample's list file and res."),
  make_option(c("-c", "--pfb"), action = "store", default = NA, type = "character",
              help = "pfb file."),
  make_option(c("-d", "--gcmodel"), action = "store", default = NA, type = "character",
              help = "gcmodel file."),
  make_option(c("-e", "--hmm"), action = "store", default = NA, type = "character",
              help = "HMM model file.")
)

## sampleID column must be changed into pure SampleID

opt = parse_args(OptionParser(option_list = option_list))

## set as parameters
## test on FA batch1
# path_dat <- "/sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/run.PennCNV/batch1/dat"
# path_main <- "/sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/run.PennCNV/batch1/res_job"

# setwd(dir = path_dat) ## no need
# system("module load penncnv") # use direct path of running command

# file_pfb <- "/sc/orga/projects/haok01a/chengh04/Food_Allergy/pipeline.callCNV/PennCNV/SNP.pfb"
# file_gcmodel <- "/sc/orga/projects/haok01a/chengh04/Food_Allergy/pipeline.callCNV/PennCNV/PCs_6_batches_12/SNP.gcmodel"
# file_hmm <- "/hpc/packages/minerva-common/penncnv/2011Jun16/lib/hhall.hmm"

path_dat     <- opt$dat
path_main    <- opt$main
file_pfb     <- opt$pfb
file_gcmodel <- opt$gcmodel
file_hmm     <- opt$hmm

if (any(is.na(c(path_dat, path_main, file_pfb, file_gcmodel, file_hmm)))) {
  stop("all parameters must be supplied. (--help for detail)")
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

  # detect_cnv.pl -test --confidence \
  # -hmm /hpc/packages/minerva-common/penncnv/2011Jun16/lib/hhall.hmm \
  # -pfb ../../SNP.pfb -gcmodel ../SNP.gcmodel -list ../list.txt \
  # -log ../FA_PCs_2_batches_3.log -out ../FA_PCs_2_batches_3.rawcnv

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









