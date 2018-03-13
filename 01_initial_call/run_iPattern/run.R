
# run-iPattern ---------------------------------------------------------------

module load R/3.0.3
module load python

export IPNBASE='/sc/orga/projects/haok01a/chengh04/shared_genomics_resources/iPattern/FA/ipn_0.581'
PYTHONPATH=$PYTHONPATH:'/sc/orga/projects/haok01a/chengh04/shared_genomics_resources/iPattern/FA/ipn_0.581/ipnlib'

/sc/orga/projects/haok01a/chengh04/shared_genomics_resources/iPattern/FA/ipn_0.581/preprocess/ilmn/ilmn_run.py \
--gender-file /sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/run.batch2/iPattern/2_1/FA_batch_2_1_gender.txt \
--bad-sample-file /sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/run.batch2/iPattern/2_1/FA_batch_2_1_bad_samples.txt \
--data-file-list /sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/run.batch2/iPattern/2_1/FA_batch_2_1_data_file.txt \
--experiment FA_batch_2_1 \
-o /sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/run.batch2/iPattern/2_1 \
--do-log --do-cleanup --noqsub

