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
- [03 CNVR](#03-CNVR)
- [04 genotype](#04-genotype)
- [05 boundary refinement](#05-boundary-refinement)
- [06 result](#06-result)
  - [compare duplicate pairs consistency rate](#compare-duplicate-pairs-consistency-rate)


## 01 initial all

prepare all BAF and LRR matrix 

call iPattern, PennCNV and QuantiSNP

### prepare chr-based matrix

Before running this script, the following data must by supplied.
1, 
```perl
perl chr_split_finalreport.pl FA_finalreport.txt \
idx_sampleID  idx_chr idx_position idx_snpName idx_LRR \
path_save_chr path_save_summary
```

### prepare data for IPQ

### call PennCNV

Here, calling PennCNV including following 5 steps:

prepare files for running PennCNV:
```sh
step.0.prepare.files.sh
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

combine all PennCNV calling results:
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
./prepare.QuantiSNP.R \
-i path/to/data/folder \
-o path/to/result/folder
```
check jobs and resubmit:
```sh
./check_QuantiSNP.R \
-d path/to/data/folder \
-r path/to/callingCNV/folder 
```

combine CNV calling results:
```sh
./combine.CNVcalling.QuantiSNP.R \
-r path/to/results/folder \
-o path/to/save/combined-result \
-n saving_name
```


### call iPattern

sample script for calling iPattern:
```sh
module load R/3.0.3
module load python

export IPNBASE='/sc/orga/projects/haok01a/chengh04/shared_genomics_resources/iPattern/FA/ipn_0.581'
PYTHONPATH=$PYTHONPATH:'/sc/orga/projects/haok01a/chengh04/shared_genomics_resources/iPattern/FA/ipn_0.581/ipnlib'

/sc/orga/projects/haok01a/chengh04/shared_genomics_resources/iPattern/FA/ipn_0.581/preprocess/ilmn/ilmn_run.py \
--gender-file /sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/run.iPattern/batch1/group1/FA_batch1_group1_gender.txt \
--bad-sample-file /sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/run.iPattern/batch1/group1/FA_batch1_group1_bad_samples.txt \
--data-file-list /sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/run.iPattern/batch1/group1/FA_batch1_group1_data_file.txt \
--experiment FA_batch1_group1 \
-o /sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/run.iPattern/batch1/group1 \
--do-log --do-cleanup --noqsub
```


## 02 batch effect

### snp-level LRR statics
 
randomly select 100000 snps

```r
randomly.select.snps.R
```

generate snps LRR matrix

```perl
perl generate.snps.LRR.matrix.pl
```
PCA

```r
PCA.R
```

### sample-level IPQ 10 statics

## 03 CNVR

contain one method and IPQ merge

## 04 genotype

genotyping for all CNVRs


## 05 boundary refinement

boundary refinement

```sh
./boundary_refinement.R -c 1 \
-r path/to/cnvr.rds -l path/to/chr-lrr-matrix \
-p path/to/snp.pfb -m path/to/chr-centromere.rds \
-g path/to/save/png -o path/to/save/detail-results \
-s path/to/rcpp
```
you need to source following script.

```r
refine_step1.cpp
```

## 06 result

summary compare results

### compare duplicate pairs consistency rate

