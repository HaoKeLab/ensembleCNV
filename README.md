# ensembleCNV

## Method Description

EnsembleCNV, which first detect CNV by aggregating the complementary strengths from multiple existing callers, followed by re-genotype and boundary refinement.  

## Table of Contents

- [01 initial all](#01-initial-all)
  - [prepare chr-based LRR matrix and BAF matrix](#prepare-chr-based-lrr-matrix-and-baf-matrix)
  - [prepare data for running IPQ](#prepare-data-for-running-IPQ)
  - [call PennCNV](#call-penncnv)
  - [call QuantiSNP](#call-quantisnp)
  - [call iPattern](#call-ipattern)
- [02 batch effect](#02-batch-effect)
  - [snp-level LRR statics](#snp-level-lrr-statics)
  - [sample-level IPQ 10 statics](#sample-level-ipq-10-statics)
- [03 create CNVR](#03-create-CNVR)
  - [create CNVR for individual CNV calling method](#create-CNVR-for-individual-CNV-calling-method)
  - [ensembleCNV](#ensenmbleCNV)
- [04 genotype](#04-genotype)
  - [split cnvrs into batches](#split-cnvr-into-batches)
  - [regenotype](#regenotype)
  - [combine prediction results](#combine-prediction-results)
- [05 boundary refinement](#05-boundary-refinement)
- [06 result](#06-result)
  - [compare duplicate pairs consistency rate](#compare-duplicate-pairs-consistency-rate)
- [test](#test)
  - [test ensemblCNV](#test-ensembleCNV)
  - [test regenotype](#test-regenotype)


## 01 initial all

prepare all BAF and LRR matrix 

call iPattern, PennCNV and QuantiSNP

### prepare chr-based matrix (LRR and BAF)

Before running this script, the following data must by supplied.
1, generate LRR and BAF (tab format) matrix from finalreport
```perl
perl finalreport_to_matrix_LRR_and_BAF.pl \
path_to_finalreport \
path_to_save_matrix_in_tab_format
```
2, tansform tab format matrix to .rds format
```sh
./tranform_from_tab_to_rds.R path_input path_output chr_start chr_end
```

### prepare data for IPQ

iPattern
```sh
perl finalreport_to_iPattern.pl \
-prefix path_to_save_ipattern_input_file \
-suffix .txt \
path_to_finalreport
```

PennCNV
```sh
perl finalreport_to_PennCNV.pl \
-prefix path_to_save_penncnv_input_file \
-suffix .txt \
path_to_finalreport
```

QuantiSNP
```sh
perl finalreport_to_QuantiSNP.pl \
-prefix path_to_save_quantisnp_input_file \
-suffix .txt \
path_to_finalreport
```

### call PennCNV

Here, calling PennCNV including following 5 steps:

prepare files containing SNP.pfb and SNp.gcmodel for running PennCNV:
```sh
step.0.prepare.files.sh contains all commands 
```

run PennCNV through submiting jobs:
```sh 
./step.1.run.PennCNV.jobs.R \
-a path/to/dat \
-b path/to/res_job \
-c path/to/SNP.pfb \
-d path/to/SNP.gcmodel \
-e path/to/penncnv/2011Jun16/lib/hhall.hmm
```

check jobs and resubmit unfinishing callings:
```sh
./step.2.check.PennCNV.R \
-a path/to/dat \
-b path/to/res_job \
-c path/to/SNP.pfb \
-d path/to/SNP.gcmodel \
-e path/to/penncnv/2011Jun16/lib/hhall.hmm
```

combine all PennCNV calling results (sample based):
```sh
perl step.3.combine.PennCNV.pl \
--in_dir path/to/res/ \
--out_dir path/to/output/
```

clean PennCNV and generate final results:
```sh
./step.4.clean.PennCNV.R \
-i path/to/result/folder \
-p path/to/SNP.pfb \
-n saving_name
```

### call QuantiSNP

Here, calling QuantiSNP including 3 steps:

prepare QuantiSNP and submit jobs:
```sh
./step.1.prepare.QuantiSNP.R \
-i path/to/data/folder \
-o path/to/result/folder
```
check jobs and resubmit:
```sh
./step.2.check_QuantiSNP.R \
-d path/to/data/folder \
-r path/to/callingCNV/folder 
```

combine CNV calling results:
running this script, you need to add "in_dir", "out_dir", "out_file" information in the script.
```sh
perl step.3.combine.QuantiSNP.pl
```


### call iPattern

sample script for calling iPattern:
```sh
script "run.R" contains all needed running command.
```


## 02 batch effect

### PCA on snp-level LRR statics from randomly select 100000 snps

```sh
./step.1.randomly.select.snp.R file_snps path_output

perl step.2.generate.snps.LRR.matrix.pl (add "file_snps_selected", "finalreport", "file_matrix_LRR")

step.3.pca.new.R ( add "filename_matrix", "path_input")
``` 

### PCA on sample-level iPattern, PennCNV and QuantiSNP generated 10 statics

```sh

generate iPattern, PennCNV and QuantiSNP calling sample level statics data using step.1.generate.data.R script

do PCA using step.2.pca.R
```

## 03 create CNVR

Here, create CNVR for both individual CNV calling method and ensembleCNV.

First, CNV calling results (.rds format) from iPattern, PennCNV and QuantiSNP 
```sh
step.1.data.R ( "path_output", "file_ipattern", "file_penncnv", "file_quantisnp" )
```

Second, create CNVR
individual method:
```sh 
./step.2.create.CNVR.IPQ.R --help for detail
```
ensembleCNV method:
```sh
./step.2.ensembleCNV.R --help for detail
```

Third, generate matrix for individual CNV calling method:
```sh
./step.3.generate.matrix.R file_cnv cnvCaller path_output
```

## 04 genotype

genotyping for all CNVRs containing two main steps:

split all cnvrs generated from ensembleCNV step into chromosome based batches.
```sh
./step.1.split.cnvrs.into.batches.R --help for detail
```

regenotype CNVRs in one batch:
(1) path_data contains:
samples_QC.rds (from PennCNV with columns: LRR_mean and LRR_sd )
duplicate.pairs.rds (two columns: sample1.name sample2.name )
SNP.pfb (from PennCNV)
cnvs_step3_clean.rds (from ensembleCNV step)
cnvrs_batch_annotated.rds (cnvrs adding columns: batch)
(2) path_matrix (LRR and BAF folder to save matrix data)
(3) path_sourcefile (all scripts need to be source )
(4) path_result (save all regenotype results)

```sh
sample code:
./step.2.regenotype.each.chr.each.batch.R \
-c 1 -b 1 -t 0 -p path_data -o path_result -m path_matrix -s path_sourcefile

./step.2.regenotype.each.chr.each.batch.R --help for detail

```

combine all sample-based regenotype results.
and, generate mat_CN.rds (matrix of regenotype copy number),
matrix_GQ.rds (matrix of regenotype gq score), 
CNVR_ID.rds (rownames of matrix),
Sample_ID.rds( columns of matrix).

explation:
path_cnvr (with cnvrs_annotated_batch.rds)
path_pred (with chr-batch-based regenotype results)
path_res (save results: matrix_CN.rds matrix_GQ.rds)
```sh
./step.5.prediction.results.R n.samples path_cnvr path_pred pred_res
```

## 05 boundary refinement

There are 5 steps in boundary refinement, as following:

All scripts are in folder 05_boundary_refinement.

The main part is script named as step.2.boundary_refinement.R:
```sh
./step.2.boundary_refinement.R --help for detail
```

## 06 result

summary compare results between all CNV calling methods with ensembleCNV method.
copy all following files to path_input:
dup_samples.rds with columns: sample1.name, sample1.name
matrix_iPattern.rds; matrix_PennCNV.rds; matrix_QuantiSNP.rds; 
matrix_IPQ_intersect.rds; matrix_IPQ_union.rds; matrix_ensembleCNV.rds

### compare duplicate pairs consistency rate, call rate based on sample and CNVR.

```sh
./compare.dups.consistency.R path_input cohort_name path_output
```
## test
here, we supply a samll test example for user to test.

### Test ensembleCNV

### Test regenotype


