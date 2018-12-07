#!/urs/bin/env Rscript

## The script was used to prepare auxiliary input files for iPattern
## the auxiliary input files will be stored at ${WKDIR}/01_initial_call/run_iPattern/data_aux

args <- commandArgs( trailingOnly = TRUE )

## working directory
path_wkdir <- args[1]

## project name for running iPattern
project_name <- args[2]

## path to the directory for running iPattern
path_run_ipattern <- file.path(path_wkdir, "01_initial_call", "run_iPattern")

##--------------------------------------------------------------------------------
## 1) data_file: list of splitted final report files for each sample
## the directory contains the input files prepared by finalreport_to_iPattern.pl
path_ipattern_prepare_data <- file.path(path_run_ipattern, "data")
fls_all <- list.files(path = path_ipattern_prepare_data, pattern = ".txt$", full.names = TRUE)

data_file <- data.frame(data_file = fls_all, stringsAsFactors = FALSE)
write.table( data_file, file = file.path( path_run_ipattern, "data_aux", 
                                          paste0(project_name, "_data_file.txt")),
             sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


##--------------------------------------------------------------------------------
## 2) gener_file: tab-delimied file which lists geneder information for each sample
## the file consists of two columns, Sample_ID and Gender,
## which may be retrieved from Samples_Table.txt (see Data section of ensembleCNV README.md)
## Samples_Table.txt is supposed to be at ${WKDIR}/data
## the gender_file does NOT have column names in the header, for example
# Sample_1	M
# Sample_2	F
# Sample_3	F

gender_file <- read.delim(file = file.path(path_wkdir, "data", "Samples_Table.txt"), as.is = TRUE)
gender_file$Gender <- toupper( substr(gender_file$Gender, 1, 1) )
write.table( gender_file, file = file.path( path_run_ipattern, "data_aux", 
                                            paste0(project_name, "_gender_file.txt")),
             sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


##--------------------------------------------------------------------------------
## 3) bad_samples: file lists sample IDs of poor quality to be excluded from iPattern analysis, for example
# bad_sample_1
# bad_sample_2
# bad_sample_3

## We prepare an empty file. The user can type in bad samples.
write.table(NULL, file = file.path( path_run_ipattern, "data_aux",
                                    paste0(project_name, "_bad_samples.txt")),
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

cat("Processing is completed.\n")
cat("Three files are generated:\n")
cat(file.path( path_run_ipattern, "data_aux", paste0(project_name, "_data_file.txt")), "\n")
cat(file.path( path_run_ipattern, "data_aux", paste0(project_name, "_gender_file.txt")), "\n")
cat(file.path( path_run_ipattern, "data_aux", paste0(project_name, "_bad_samples.txt")), "\n")

