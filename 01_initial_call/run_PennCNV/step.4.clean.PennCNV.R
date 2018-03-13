#!/hpc/packages/minerva-common/R/3.3.1/lib64/R/bin/Rscript --vanilla

suppressPackageStartupMessages(require(optparse))

option_list <- list(
  make_option(c("-i", "--input"), action = "store", default = NA, type = "character",
              help = "input path for combined PennCNV result."),
  make_option(c("-p", "--pfb"), action = "store", default = NA, type = "character",
              help = "pfb file"),
  make_option(c("-n", "--name"), action = "store", default = NA, type = "character",
              help = "rawcnv filename")
)

opt <- parse_args(OptionParser(option_list = option_list))

path_input  <- opt$input
file_pfb    <- opt$pfb
name_rawcnv <- opt$name

if (any(is.na(c(path_input, file_pfb, name_rawcnv)))) {
  stop("all parameters must be supplied.( --help for detail )")
}

# path_res <- "/sc/orga/projects/haok01a/chengh04/Food_Allergy/run.PennCNV.jobs/res"
# path_dat <- "/sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/run.PennCNV/batch1/dat"
# 
# sample_files <- list.files(path = path_dat) ## all dat files  ## form raw data files
# samples  <- gsub(pattern = ".txt$", replacement = "", x = sample_files)
# n.sample <- length(samples)

# clean CNV ---------------------------------------------------------------

setwd(dir = path_input)
path_clean <- path_input
name_project <- name_rawcnv

file_rawcnv <- paste(name_project, "rawcnv", sep = ".")
## make the default pfb file
# file_pfb <- "/sc/orga/projects/haok01a/chengh04/Food_Allergy/pipeline.callCNV/PennCNV/SNP.pfb"
file_pfb <- file_pfb

n_rawcnv <- system(paste("cat", file_rawcnv, "| wc -l"), intern = TRUE)
n_rawcnv <- as.integer(n_rawcnv)

cat("CNV number before clean:", n_rawcnv, "\n")

flag = 0
idx = 1
cnv1_in <- file_rawcnv
while(flag == 0) {
  
  n_rawcnv <- as.integer(system(paste("cat", cnv1_in, "|", "wc -l"), intern = TRUE))
  cnv1_out <- paste(name_project, idx, "rawcnv", sep = ".")
  
  cmd1 <- paste("/hpc/packages/minerva-common/penncnv/2011Jun16/bin/clean_cnv.pl",
                "combineseg", cnv1_in, "--signalfile", file_pfb, 
                "--fraction 0.2", "--bp >", cnv1_out)
  
  cat("Start run Time:", idx, cmd1, "...\n")
  system(cmd1)
  cat("End run ......\n")
  
  cmd2 <- paste("cat", cnv1_out, "|", "wc -l")
  n_newcnv <- system(cmd2, intern = TRUE)
  n_newcnv <- as.integer(n_newcnv)
  
  cat("raw number:", n_rawcnv, "\n")
  cat("new number:", n_newcnv, "\n")
  
  if (n_rawcnv == n_newcnv) {
    flag = 1
  } else {
    cnv1_in <- cnv1_out
    idx <- idx + 1
  }
  
}

## transform format of CNV results
# convert_cnv.pl --intype penncnv --outtype tab Valentina_112samples.2.rawcnv > Valentina_112samples.txt
cnv_penncnv <- paste(name_project, idx, "rawcnv", sep = ".")
cnv_tab <- paste(name_project, "txt", sep = ".")
cat("transform format of CNV final result.\n")
cmd.transform <- paste("/hpc/packages/minerva-common/penncnv/2011Jun16/bin/convert_cnv.pl",
                       "--intype", "penncnv", "--outtype", "tab", cnv_penncnv, ">", cnv_tab)
system(cmd.transform)

## extract individual level statistics for QC 
cat("extract individual level statistics for QC.\n")
cnv_log <- paste(name_project, "log", sep = ".")
cnv_qc <- paste0(name_project, "_qc.txt")
cmd.extract <- paste("/hpc/packages/minerva-common/penncnv/2011Jun16/bin/filter_cnv.pl", cnv_penncnv,
                     "-qclogfile", cnv_log, "-qcsumout", cnv_qc, ">", "step4.txt")
system(cmd.extract)

# Change SampleID column information --------------------------------------
# path to pure SampleID
## CNV
dat_CNV <- read.table(file = cnv_tab, sep = "\t",
                      header = FALSE, comment.char = "", check.names = FALSE, stringsAsFactors = FALSE)
samples_path <- dat_CNV$V5
sampleIDs <- unlist(lapply(1:length(samples_path), FUN = function(k) {
  sample1 <- samples_path[k]
  str1 <- unlist(strsplit(sample1, split = "/", fixed = TRUE))
  str1[length(str1)]
}))
dat_CNV$V5 <- sampleIDs  ## change

write.table(dat_CNV, file = paste0(name_project, "_new.txt"),
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

## Sample_Stat
dat_Sample_Stat <- read.table(file = cnv_qc, sep = "\t",
                              header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
files <- dat_Sample_Stat$File
files_new <- unlist(lapply(1:length(files), FUN = function(k) {
  file1 <- files[k]
  str1  <- unlist(strsplit(file1, split = "/", fixed = TRUE))
  str1[length(str1)]
}))

dat_Sample_Stat$File <- files_new

write.table(dat_Sample_Stat, file = paste0(name_project, "_qc_new.txt"),
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)







