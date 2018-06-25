


path_run_ipattern <- ""
# perpare -----------------------------------------------------------------

nm_prefix = ""

## data_file
path_ipattern_prepare_data <- ""
fls_all <- list.files(path = path_ipattern_prepare_data, pattern = ".txt$")

data_file <- data.frame(data_file = fls_all, stringsAsFactors = FALSE)
write.table( data_file, file = file.path( path_run_ipattern, paste(nm_prefix, "data_file.txt", sep = "_")),
             sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

## geneder
## without columnName: Sample_ID and gender
# EFZ160736A      M
# EFZ80669A       F
# EFZ157712A      F

## bad_samples
# if there are no bad samples to exclude when running ipattern
write.table(NULL, file = file.path( path_run_ipattern, paste(nm_prefix, "bad_samples.txt", sep = "_")),
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# run-iPattern ---------------------------------------------------------------

module load R/3.0.3 ## other will cause wrong
module load python

path_to_ipattern=""
export IPNBASE='/sc/orga/projects/haok01a/chengh04/shared_genomics_resources/iPattern/FA/ipn_0.581'
PYTHONPATH=$PYTHONPATH:'/sc/orga/projects/haok01a/chengh04/shared_genomics_resources/iPattern/FA/ipn_0.581/ipnlib'

INPUT_PATH=""

${path_to_ipattern}/ipn_0.581/preprocess/ilmn/ilmn_run.py \
--gender-file ${INPUT_PATH}FA_batch_2_1_gender.txt \
--bad-sample-file ${INPUT_PATH}FA_batch_2_1_bad_samples.txt \
--data-file-list ${INPUT_PATH}FA_batch_2_1_data_file.txt \
--experiment FA_batch_2_1 \
-o /sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/run.batch2/iPattern/2_1 \
--do-log --do-cleanup --noqsub

