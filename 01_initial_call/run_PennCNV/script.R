

## run PennCNV FA 3 batches ------------------------------------------------------------

# batch1
./step.1.run.PennCNV.jobs.R \
-a /sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/run.PennCNV/batch1/dat \
-b /sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/run.PennCNV/batch1/res_job \
-c /sc/orga/projects/haok01a/chengh04/Food_Allergy/pipeline.callCNV/PennCNV/SNP.pfb \
-d /sc/orga/projects/haok01a/chengh04/Food_Allergy/pipeline.callCNV/PennCNV/PCs_6_batches_12/SNP.gcmodel \
-e /hpc/packages/minerva-common/penncnv/2011Jun16/lib/hhall.hmm

# batch2
./step.1.run.PennCNV.jobs.R \
-a /sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/run.PennCNV/batch2/dat \
-b /sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/run.PennCNV/batch2/res_job \
-c /sc/orga/projects/haok01a/chengh04/Food_Allergy/pipeline.callCNV/PennCNV/SNP.pfb \
-d /sc/orga/projects/haok01a/chengh04/Food_Allergy/pipeline.callCNV/PennCNV/PCs_6_batches_12/SNP.gcmodel \
-e /hpc/packages/minerva-common/penncnv/2011Jun16/lib/hhall.hmm

# batch3
./step.1.run.PennCNV.jobs.R \
-a /sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/run.PennCNV/batch3/dat \
-b /sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/run.PennCNV/batch3/res_job \
-c /sc/orga/projects/haok01a/chengh04/Food_Allergy/pipeline.callCNV/PennCNV/SNP.pfb \
-d /sc/orga/projects/haok01a/chengh04/Food_Allergy/pipeline.callCNV/PennCNV/PCs_6_batches_12/SNP.gcmodel \
-e /hpc/packages/minerva-common/penncnv/2011Jun16/lib/hhall.hmm

## check PennCNV results for each batch -----------------------------------------------
# batch1
./step.2.check.PennCNV.R \
-a /sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/run.PennCNV/batch1/dat \
-b /sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/run.PennCNV/batch1/res_job \
-c /sc/orga/projects/haok01a/chengh04/Food_Allergy/pipeline.callCNV/PennCNV/SNP.pfb \
-d /sc/orga/projects/haok01a/chengh04/Food_Allergy/pipeline.callCNV/PennCNV/PCs_6_batches_12/SNP.gcmodel \
-e /hpc/packages/minerva-common/penncnv/2011Jun16/lib/hhall.hmm

# batch2
./step.2.check.PennCNV.R \
-a /sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/run.PennCNV/batch2/dat \
-b /sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/run.PennCNV/batch2/res_job \
-c /sc/orga/projects/haok01a/chengh04/Food_Allergy/pipeline.callCNV/PennCNV/SNP.pfb \
-d /sc/orga/projects/haok01a/chengh04/Food_Allergy/pipeline.callCNV/PennCNV/PCs_6_batches_12/SNP.gcmodel \
-e /hpc/packages/minerva-common/penncnv/2011Jun16/lib/hhall.hmm

# batch3
./step.2.check.PennCNV.R \
-a /sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/run.PennCNV/batch3/dat \
-b /sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/run.PennCNV/batch3/res_job \
-c /sc/orga/projects/haok01a/chengh04/Food_Allergy/pipeline.callCNV/PennCNV/SNP.pfb \
-d /sc/orga/projects/haok01a/chengh04/Food_Allergy/pipeline.callCNV/PennCNV/PCs_6_batches_12/SNP.gcmodel \
-e /hpc/packages/minerva-common/penncnv/2011Jun16/lib/hhall.hmm

## combine PennCNV results -----------------------------------------------
# batch1
perl step.3.combine.PennCNV.pl \
--in_dir /sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/run.PennCNV/batch1/res_job/res/ \
--out_dir /sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/run.PennCNV/batch1/res_job/

# batch2
perl step.3.combine.PennCNV.pl \
--in_dir /sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/run.PennCNV/batch2/res_job/res/ \
--out_dir /sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/run.PennCNV/batch2/res_job/

# batch3
perl step.3.combine.PennCNV.pl \
--in_dir /sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/run.PennCNV/batch3/res_job/res/ \
--out_dir /sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/run.PennCNV/batch3/res_job/

## clean PennCNV and generate CNV and Stat ------------------------------------

## batch1
./step.4.clean.PennCNV.R \
-i /sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/run.PennCNV/batch1/res_job \
-p /sc/orga/projects/haok01a/chengh04/Food_Allergy/pipeline.callCNV/PennCNV/SNP.pfb \
-n CNV.PennCNV

## batch2
./step.4.clean.PennCNV.R \
-i /sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/run.PennCNV/batch2/res_job \
-p /sc/orga/projects/haok01a/chengh04/Food_Allergy/pipeline.callCNV/PennCNV/SNP.pfb \
-n CNV.PennCNV

## batch3
./step.4.clean.PennCNV.R \
-i /sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/run.PennCNV/batch3/res_job \
-p /sc/orga/projects/haok01a/chengh04/Food_Allergy/pipeline.callCNV/PennCNV/SNP.pfb \
-n CNV.PennCNV



## other not jobs results 
./step.4.clean.PennCNV.R \
-i /sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/run.PennCNV/batch1 \
-p /sc/orga/projects/haok01a/chengh04/Food_Allergy/pipeline.callCNV/PennCNV/SNP.pfb \
-n FA_batch1







