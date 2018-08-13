## The script was used to run iPattern on Minerva high performance cluster.
## You need to modifiy it according to the system you are using if you would like to use it.

## working directory
path_run_ipattern <- ""

# perpare input data-------------------------------------------------------
## project specific prefix
nm_prefix = ""

## 1) data_file: list of splitted final report files for each sample
## the directory contains the input files prepared by finalreport_to_iPattern.pl
path_ipattern_prepare_data <- ""
fls_all <- list.files(path = path_ipattern_prepare_data, pattern = ".txt$")

data_file <- data.frame(data_file = fls_all, stringsAsFactors = FALSE)
write.table( data_file, file = file.path( path_run_ipattern, paste(nm_prefix, "data_file.txt", sep = "_")),
             sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

## 2) gener_file: tab-delimied file which lists geneder information for each sample
## the file consists of two columns, Sample_ID and gender,
## withou column names in the header, for example
# Sample_1	M
# Sample_2	F
# Sample_3	F

## 3) bad_samples: file lists sample IDs of poor quality to be excluded from iPattern analysis, for exampl
# bad_sample_1
# bad_sample_2
# bad_sample_3

## if there are no bad samples to be excluded, just prepare an empty file
write.table(NULL, file = file.path( path_run_ipattern, paste(nm_prefix, "bad_samples.txt", sep = "_")),
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


# run-iPattern on linux command line ------------------------------------------

path_to_ipattern="" ## where iPattern is installed
export IPNBASE="$path_to_ipattern/ipn_0.581"
PYTHONPATH=$PYTHONPATH:"$path_to_ipattern/ipn_0.581/ipnlib"

INPUT_PATH="" ## where gender_file.txt, bad_file.txt and data_file.txt are located

${path_to_ipattern}/ipn_0.581/preprocess/ilmn/ilmn_run.py \
--gender-file ${INPUT_PATH}/gender_file.txt \
--bad-sample-file ${INPUT_PATH}/bad_samples.txt \
--data-file-list ${INPUT_PATH}/data_file.txt \
--experiment nm_prefix \
-o path_run_ipattern \
--do-log --do-cleanup --noqsub

