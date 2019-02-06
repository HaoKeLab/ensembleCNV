# ensembleCNV

## Method Description

EsembleCNV is a novel ensemble learning framework to detect and genotype copy number variations (CNVs) using single nucleotide polymorphism (SNP) array data. EnsembleCNV a) identifies and eliminates batch effects at raw data level; b) assembles individual CNV calls into CNV regions (CNVRs) from multiple existing callers with complementary strengths by a heuristic algorithm; c) re-genotypes each CNVR with local likelihood model adjusted by global information across multiple CNVRs; d) refines CNVR boundaries by local correlation structure in copy number intensities; e) provides direct CNV genotyping accompanied with confidence score, directly accessible for downstream quality control and association analysis. 

More details can be found in the manuscript:

Zhongyang Zhang, Haoxiang Cheng, Xiumei Hong, Antonio F. Di Narzo, Oscar Franzen, Shouneng Peng, Arno Ruusalepp, Jason C. Kovacic, Johan LM Bjorkegren, Xiaobin Wang, Ke Hao (2019) EnsembleCNV: An ensemble machine learning algorithm to identify and genotype copy number variation using SNP array data. Nucleic Acids Research, doi: https://doi.org/10.1093/nar/gkz068

The scripts are in its beta version. Please report bugs and issues or provide suggestions at [here](https://github.com/HaoKeLab/ensembleCNV/issues).

## Table of Contents

- [Installation](#installation)
- [Data](#data)
- [1 Initial call](#1-initial-call)
  - [Prepare chromosome-wise LRR and BAF matrices for CNV genotyping](#prepare-chromosome-wise-lrr-and-baf-matrices-for-cnv-genotyping)
  - [Prepare data for individual CNV callers](#prepare-data-for-individual-cnv-callers)
- [2 Batch effect](#2-batch-effect)
  - [PCA on raw LRR data](#pca-on-raw-lrr-data)
  - [PCA on summary statistics](#pca-on-summary-statistics)
- [3 Create CNVR](#3-create-cnvr)
- [4 CNV genotyping for each CNVR](#4-cnv-genotyping-for-each-cnvr)
- [5 Boundary refinement](#5-boundary-refinement)
- [6 Performance assessment](#6-performance-assessment)
- [Example](#example)
  - [Create CNVR](https://github.com/HaoKeLab/ensembleCNV/tree/master/example/example_create_CNVR)
  - [CNV genotyping](https://github.com/HaoKeLab/ensembleCNV/tree/master/example/example_CNV_genotype)
  - [Boundary refinement](https://github.com/HaoKeLab/ensembleCNV/tree/master/example/example_boundary_refinement)

## Installation

### Requirements

- R (3.3.1+) (https://www.r-project.org/) with the following packages:
  - cowplot (0.9.2+)
  - data.table (1.10.4-3+)
  - dplyr (0.7.4+)
  - ggplot2 (3.0.0+)
  - gridExtra (2.3+)
  - mclust (5.4+)
  - mixtools (1.1.0+)
  - modeest (2.1+)
  - optparse (1.3.2+)
  - pheatmap (1.0.8+)
  - plyr (1.8.4+)
  - RColorBrewer (1.1-2+) 
  - Rcpp (0.12.17+)
  - RcppArmadillo (0.7.500.0.0+)
  - tibble (1.4.2+)
- Perl (5.10.1+) (https://www.perl.org/)

Note: 

- The scripts of ensembleCNV have been developed and tested on [Minerva](https://hpc.mssm.edu/), a high-performance multi-node linux cluster (CentOS 6.9) wtih LSF (Load Sharing Facility). Part of the scripts used for job submission to parallelize the computation needs to be adjusted to your specific computational environment (see below sections [CNV genotyping for each CNVR](#4-cnv-genotyping-for-each-cnvr) and [Boundary refinement](#5-boundary-refinement)). Running the whole pipeline on a Linux laptop/desktop may be possible, but may take much longer time and sometimes may be computationally prohibitive, especially for large projects with thousands to tens of thousands of samples (i.e., typical GWAS data).

- Please be advised that ensembleCNV is designed to detect and genotype CNVs on a relatively large cohort usually consisting of at least a few hundred samples. In particular, the steps [creating CNVR](#3-create-cnvr), [CNV genotyping](#4-cnv-genotyping-for-each-cnvr), and [boundary refinement](#5-boundary-refinement) require relatively large sample size to achieve a reasonable reproducibility and accuracy. Results generated from only a few samples are not valid.  

- Each of the 3 third-source CNV callers has its own installation requirements. Please refer to their specific installation instructions: [iPattern](https://github.com/HaoKeLab/ensembleCNV/tree/master/01_initial_call/run_iPattern), [PennCNV](https://github.com/HaoKeLab/ensembleCNV/tree/master/01_initial_call/run_PennCNV) and [QuantiSNP](https://github.com/HaoKeLab/ensembleCNV/tree/master/01_initial_call/run_QuantiSNP).

### Installation

```sh
git clone https://github.com/HaoKeLab/ensembleCNV
```
The scripts in the installation folder `ensembleCNV` are organized step by step, with scripts for each step located in each individual subfolder. For a new project, we recommend the user make a copy of the original installation folder `ensembleCNV` in the working directory and keep the folder structure to organize the data and analysis workflow. We prepared a shell script for creating new project.
```sh
ENSEMBLECNV=</path/to/ensembleCNV>
cd $ENSEMBLECNV
chmod +x create_new_project.sh

WKDIR=</path/to/working_directory>
./create_new_project.sh $WKDIR
```

Please go through the detailed step-by-step instructions as follows when using ensembleCNV for the first time. For your reference, we also prepared a template script [run_ensembleCNV_template.sh](https://github.com/HaoKeLab/ensembleCNV/blob/master/run_ensembleCNV_template.sh) to facilitate running the whole pipeline.

## Data

### Final report

The raw data comes from the [final report](http://jp.support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/genomestudio/genomestudio-2-0/genomestudio-genotyping-module-v2-user-guide-11319113-01.pdf) generated by Illumina [GenomeStudio](https://support.illumina.com/array/array_software/genomestudio.html). We recommend the users follow the [protocol](https://www.ncbi.nlm.nih.gov/pubmed/25321409) to process the raw Illumina SNP array data with GenomeStudio. With the GenomeStudio, the exported final report text file is supposed to include at minimum the following 10 columns:
  - Sample ID
  - SNP Name 
  - Chr
  - Position
  - Allele 1 - Forward (or Allele 1 - Top) (used by iPattern)
  - Allele 2 - Forward (or Allele 2 - Top) (used by iPattern)
  - X (used by iPattern)
  - Y (used by iPattern)
  - Log R Ratio (used by PennCNV, QuantiSNP, and ensembleCNV)
  - B Allele Freq (used by PennCNV, QuantiSNP, and ensembleCNV)

### Sample table

The users need to prepare a project-specific tab-delimited sample table with `Sample_ID` and `Gender` information for each sample. Note: Please name the table header exactly as `Sample_ID` and `Gender`. For example,
```
Sample_ID Gender
sample1 Female
sample2 Male
sample3 Male
...
```
The gender information is required by QuanitSNP and iPattern rather than ensembleCNV. Such table may have been already prepared by the investigators (i.e., this is typically the case for GWAS). Another option is to export the sample table with GenomeStudio, which has a build-in function to estimate gender if gender information is not provided by the investigators. Please refer to [GenomeStudio manual](http://jp.support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/genomestudio/genomestudio-2-0/genomestudio-genotyping-module-v2-user-guide-11319113-01.pdf) for details.

### Centromere position

We have prepared a tab-delimited table [centromere_hg19.txt](https://github.com/HaoKeLab/ensembleCNV/blob/master/example/example_create_CNVR/data/centromere_hg19.txt) for the centromere position (hg19) of each chromosome. Centromere positions for other assemblies can be extracted from corresponding Chromosome Band tables from UCSC genome browser at [here](https://genome.ucsc.edu/cgi-bin/hgTables).

### Duplicate pairs [optional]

If duplicated samples (either technical duplicates or monozygotic twins) are available, we can use the information for quality control and deciding the GQ score threshold in the step [performance assessment](#6-performance-assessment). Please refer to the [manuscript](https://doi.org/10.1093/nar/gkz068) for details. We provdie an example table [duplicate_pairs.txt](https://github.com/HaoKeLab/ensembleCNV/blob/master/example/example_CNV_genotype/data/duplicate_pairs.txt). Note: Please name the table header exactly as `sample1.name` and `sample2.name`.

Before running the analysis, please put (or create symbolic link to) the final report file (e.g., named "final_report.txt"), sample table (exactly named as "Samples_Table.txt"), centromere position table (e.g., named "centromere_hg19.txt"), and duplicates table [optional] (exactly named as "duplicate_pairs.txt") in the folder `${WKDIR}/data`.

## 1 Initial call

The pipeline begins with running inividual CNV callers, including [iPattern](https://www.ncbi.nlm.nih.gov/pubmed/?term=21552272), [PennCNV](http://penncnv.openbioinformatics.org/en/latest/), and [QuantiSNP](https://sites.google.com/site/quantisnp/), to make initial CNV calls. Before that, the final report data needs to be converted into proper format required by ensembleCNV as well as inividual CNV callers.

### Prepare chromosome-wise LRR and BAF matrices for CNV genotyping

We provide [perl scripts](https://github.com/HaoKeLab/ensembleCNV/tree/master/01_initial_call/finalreport_to_matrix_LRR_and_BAF) to extract LRR and BAF information from final report, combine them across individuals and divide them by chromsomes. 


(1) Create LRR and BAF (tab delimited) matrices from final report
```sh
perl ${WKDIR}/01_initial_call/finalreport_to_matrix_LRR_and_BAF/finalreport_to_matrix_LRR_and_BAF.pl \
${WKDIR}/data/final_report.txt \
${WKDIR}/01_initial_call/finalreport_to_matrix_LRR_and_BAF
```
Note: 
- This perl script is designed to process a final report file with multiple samples. Although the code can accommodate a final report file with one or a few samples, the output files in the case are not meaningful for downstream analysis. The script will issue a warning if sample size is too small for subsequent steps.

- In addtion to creating LRR and BAF matrices, the perl script also screens the whole final report file and performs sanity check. It will stop and report errors if there is any discripancy in the number of probes between samples. It will issue warnings if any columns required by ensembleCNV or the three CNV callers are missing.

(2) Tansform tab-delimited text file to .rds format for quick loading in R
```sh
Rscript ${WKDIR}/01_initial_call/finalreport_to_matrix_LRR_and_BAF/transform_from_tab_to_rds.R \
--input ${WKDIR}/01_initial_call/finalreport_to_matrix_LRR_and_BAF \
--output ${WKDIR}/01_initial_call/finalreport_to_matrix_LRR_and_BAF/RDS \
--startChr <INT> \ ## default: 1
--endChr <INT>  ## default: 22
```
The parameters `--startChr` and `--endChr` indicate the range of chromosomes (1 <= startChr <= endChr <= 22) to be processed. When `--startChr` and `--endChr` are not specified, all the autosomal chromosomes (i.e., Chr 1 ~ 22) will be processed by default. If you are interested in CNVs in a particular chromosome, e.g., chr 3, set `--startChr 3 --endChr 3`.

Note: In current version, we focus on detecting and genotyping CNVs in autosomal chromosomes. Functionality for analyzing CNVs of sex chromosomes is part of our future development plans.

When finishing running the scripts, there will be two folders `LRR` and `BAF` created under the folder specified by `--output`. In `LRR` (`BAF`) folder, you will see LRR (BAF) matrices stored in `matrix_chr_*_LRR.rds` (`matrix_chr_*_BAF.rds`) for each chromosome respectively. In the matrix, each row corresponds to a sample while each column a SNP. The data will be later used for CNV genotyping for each CNVR.

In addition, a text file named "SNP_pos.txt" with `Name`, `Chr`, and `Position` information of each probe will be generated at `${WKDIR}/01_initial_call/finalreport_to_matrix_LRR_and_BAF` and used in downstream analysis.

### Prepare data for individual CNV callers

We provide [perl scripts](https://github.com/HaoKeLab/ensembleCNV/tree/master/01_initial_call/prepare_IPQ_input_file) to extract information from final report and convert the data into the formatted input files required by iPattern, PennCNV and QuantiSNP. The required columns in final report will be retrived for each CNV caller and split into individual tab-delimited text files, each for one sample. 

#### iPattern
```sh
perl ${WKDIR}/01_initial_call/prepare_IPQ_input_file/finalreport_to_iPattern.pl \
-prefix ${WKDIR}/01_initial_call/run_iPattern/data/ \
-suffix .txt \
${WKDIR}/data/final_report.txt
```

#### PennCNV
```sh
perl ${WKDIR}/01_initial_call/prepare_IPQ_input_file/finalreport_to_PennCNV.pl \
-prefix ${WKDIR}/01_initial_call/run_PennCNV/data/ \
-suffix .txt \
${WKDIR}/data/final_report.txt
```

#### QuantiSNP
```sh
perl ${WKDIR}/01_initial_call/prepare_IPQ_input_file/finalreport_to_QuantiSNP.pl \
-prefix ${WKDIR}/01_initial_call/run_QuantiSNP/data/ \
-suffix .txt \
${WKDIR}/data/final_report.txt
```

To run each individual CNV caller, we provide complimentary scripts for [iPattern](https://github.com/HaoKeLab/ensembleCNV/tree/master/01_initial_call/run_iPattern), [PennCNV](https://github.com/HaoKeLab/ensembleCNV/tree/master/01_initial_call/run_PennCNV) and [QuantiSNP](https://github.com/HaoKeLab/ensembleCNV/tree/master/01_initial_call/run_QuantiSNP). We encourage users to consult with the original documents of these methods for more details. 


## 2 Batch effect

Two orthogonal signals can be used to identify batch effects in CNV calling: (i) Batch effects may be reflected in the first two or three PCs when principle component analysis (PCA) is applied on the raw LRR matrix. We randomly select 100,000 probes and apply PCA to the down-sampled matrix to save computational time. (ii) Along with CNV calls, the three detection methods generate sample-wise summary statistics, such as standard deviations (SD) of LRR, SD of BAF, wave factor in LRR, BAF drift, and the number of CNVs detected, reflecting the quality of CNV calls at the sample level.  Since these quantities are highly correlated among themselves and between methods, we also use PCA to summarize their information.  By examining the first two or three PCs visualized in scatter plots, we can identify sample outliers or batches that deviate from the majority of the normally behaved samples. 

Note: While isolated outliers should be excluded from downstream analysis, if batch effects are identified, the users need to re-normalize the samples within each outstanding batch with Genome Studio respectively. The initial CNV calling step with individual CNV callers need to be performed again on the updated data. The re-called CNVs will be combined with the remaining call set of good quality. Also, the chromosome-wise LRR and BAF matrices prepared for CNV genotyping need to be updated as done in the [initial step](#prepare-chromosome-wise-lrr-and-baf-matrices-for-cnv-genotyping).

### PCA on raw LRR data

This analysis is implemented in the following 3 steps.

(1) Randomly select 100,000 SNPs based on information from "SNP_pos.txt".
```sh
Rscript ${WKDIR}/02_batch_effect/PCA_on_LRR/step.1.down.sampling.R \
${WKDIR}/01_initial_call/finalreport_to_matrix_LRR_and_BAF/SNP_pos.txt \
${WKDIR}/02_batch_effect/PCA_on_LRR    ## path to snps.down.sample.txt
```

(2) Extract LRR values at the list of randomly selected SNPs across individuals from final report file generated by Genome Studio.
```sh
perl ${WKDIR}/02_batch_effect/PCA_on_LRR/step.2.LRR.matrix.pl \
${WKDIR}/02_batch_effect/PCA_on_LRR/snps.down.sample.txt \   ## the list of SNPs generated in step (1)
${WKDIR}/data/final_report.txt \                             ## generated by Genome Studio
${WKDIR}/02_batch_effect/PCA_on_LRR/LRR_matrix_for_PCA.txt   ## path to LRR_matrix_for_PCA.txt
```

(3) PCA on LRR matrix.
```sh
Rscript ${WKDIR}/02_batch_effect/PCA_on_LRR/step.3.LRR.PCA.R \
${WKDIR}/02_batch_effect/PCA_on_LRR/ \                       ## path to PCA results
${WKDIR}/02_batch_effect/PCA_on_LRR/LRR_matrix_for_PCA.txt   ## the LRR matrix generated in step (2)
``` 
When the analysis is finished, in the working directory, the first three PCs of all samples will be saved in tab-delimited text file (`LRR_PCA_res.txt`), as well as scatter plots of the first three PCs (`LRR_PCA_plots.png`). 

### PCA on summary statistics

Besides CNV calls, iPattern, PennCNV and QuantiSNP also generate 10 sample-level statistics: (a) SD of normalzied total intensity, and b) number of CNVs detected per sample from iPattern; (c) SD of LRR, (d) SD of BAF, (e) wave factor in LRR, (f) BAF drift, and (g) number of CNVs detected per sample from PennCNV; (h) SD of LRR, (i) SD of BAF, and (j) number of CNVs detected per sample from QuantiSNP. PCA can be performed in the follwoing 2 steps.

(1) Generate iPattern, PennCNV and QuantiSNP sample-level summary statistics.
```sh
Rscript ${WKDIR}/02_batch_effect/PCA_on_summary_stats/step.1.prepare.stats.R \
${WKDIR}/01_initial_call/run_iPattern/results \
${WKDIR}/01_initial_call/run_PennCNV/results \
${WKDIR}/01_initial_call/run_QuantiSNP/results \
${WKDIR}/02_batch_effect/PCA_on_summary_stats   ## summary statistics IPQ.stats.txt from iPattern, PennCNV and QuantiSNP results
```

(2) PCA on sample-level summary statistics.
```sh
Rscript ${WKDIR}/02_batch_effect/PCA_on_summary_stats/step.2.stats.PCA.R \
${WKDIR}/02_batch_effect/PCA_on_summary_stats   ## path to IPQ.stats.txt generated in step (1)
```
When the analysis is finished, in the working directory, the PCs of all samples will be saved in tab-delimited text file (`IPQ_stats_PCA_res.txt`), as well as scatter plots of the first three PCs (`IPQ_stats_PCA_plots.png`). 


## 3 Create CNVR

We define copy number variable region (CNVR) as the region in which CNVs called from different individuals by different callers substantially overlap with each other. We model the CNVR construction problem as identification of cliques (a sub-network in which every pair of nodes is connected) in a network context, where (i) CNVs detected for each individual from a method are considered as nodes; (ii) two nodes are connected when the reciprocal overlap between their corresponding CNV segments is greater than a pre-specified threshold (e.g. 30%); (iii) a clique corresponds to a CNVR in the sense that, for each CNV (node) belonging to the CNVR (clique), its average overlap with all the other CNVs of this CNVR is above a pre-specified threshold (e.g. 30%). The computational complexity for clique identification can be dramatically reduced in this special case, since the CNVs can be sorted by their genomic locations and the whole network can be partitioned by chromosome arms – CNVs from different arms never belong to the same CNVR. More details can be found in the [manuscript](https://doi.org/10.1093/nar/gkz068).

The algorithm is implemented in the following two steps.

(1) Extract CNV information from individual calls made by iPattern, PennCNV and QuantiSNP
```sh
Rscript ${WKDIR}/03_create_CNVR/step.1.CNV.data.R \
${WKDIR}/03_create_CNVR \
${WKDIR}/01_initial_call/run_iPattern/results/<project_name>_all_calls.txt \  ## <project_name> used in iPattern analysis
${WKDIR}/01_initial_call/run_PennCNV/results/CNV.PennCNV_new.txt \
${WKDIR}/01_initial_call/run_QuantiSNP/results/quantisnp.cnv \
${WKDIR}/data/Samples_Table.txt
```
After finishing this step, three tab-delimited tables for each respective method, `cnv.ipattern.txt`, `cnv.penncnv.txt`, and `cnv.quantisnp.txt`, will be generated with such fields as `Sample_ID`, `chr`, `posStart`, `posEnd`, `CNV_type`, etc. These files will be used as input in the following step (2).

(2) Merge CNV calls from individual methods into CNVRs
```sh
Rscript ${WKDIR}/03_create_CNVR/step.2.create.CNVR.R \
--icnv ${WKDIR}/03_create_CNVR/cnv.ipattern.txt \
--pcnv ${WKDIR}/03_create_CNVR/cnv.penncnv.txt \
--qcnv ${WKDIR}/03_create_CNVR/cnv.quantisnp.txt \
--snp ${WKDIR}/01_initial_call/finalreport_to_matrix_LRR_and_BAF/SNP_pos.txt \
--centromere ${WKDIR}/data/centromere_hg19.txt   ## for other assemblies, check UCSC genome browser (see above)
```
Two tab-delimited tables will be generated in this step: i) `cnvr_clean.txt` with the information for each constructed CNVR; ii) `cnv_clean.txt` with the information for each CNV calls from individual methods, including which CNVR each CNV belongs to.

We provide an [example](https://github.com/HaoKeLab/ensembleCNV/tree/master/example/example_create_CNVR) of this step corresponding to one example CNVR.

## 4 CNV genotyping for each CNVR

The initial CNV calls within a CNVR may be mixed with false positives and false negatives from the initial call set. Moreover, the baseline LRR value corresponding to normal CN status may substantially deviate from 0, violating the essential model assumptions for individual-wise CNV callers (e.g., PennCNV and QuantiSNP). To address these issues, we re-genotyped CN status per individual at each CNVR by a locally fitted likelihood model, with information from other CNVRs borrowed for the initialization of model parameters. Both the LRR and BAF signals from SNP probes and the LRR signal from CNV probes within a particular CNVR were used for model fitting. More details can be found in the [manuscript](https://doi.org/10.1093/nar/gkz068).

In current implementation, `CNV.genotype.one.chr.one.batch.R` is the main script that performs CNV genotyping on one batch of CNVRs at a time. It loads the R functions in the subdirectory [scripts](https://github.com/HaoKeLab/ensembleCNV/tree/master/04_CNV_genotype/scripts). The main script can genotype one CNVR at a time (i.e., a batch includes only one CNVR). Please see this [example](https://github.com/HaoKeLab/ensembleCNV/tree/master/example/example_CNV_genotype). 

Note: 
- We split the core scripts and the platform-specific workflow scripts (see below) in order to provide flexibity for users to develop their own workflow scripts specific to their particular platform, especially when a high-performance cluster is not accessible to the users. 

- We provide workflow scripts for a cluster environment, where CNVRs within different chromosomes are processed in parallel, and CNVRs within the same chromosomes are further grouped into batches for additional level of parallelization. Relevant R scripts can be found [here](https://github.com/HaoKeLab/ensembleCNV/tree/master/04_CNV_genotype). 

Before running the script below, the following files generated in previous steps need to be copied (or linked) into the `${WKDIR}/04_CNV_genotype/data` directory, and named exactly as follows:

  - `SNP.pfb -> ${WKDIR}/01_initial_call/run_PennCNV/data_aux/SNP.pfb` (prepared when running PennCNV; containing the column of PFB (Population Frequency of B allele) used in modeling the likelihood of BAF data)
  - `cnvr_clean.txt -> ${WKDIR}/03_create_CNVR/cnvr_clean.txt` (generated in "create CNVR" step)
  - `cnv_clean.txt -> ${WKDIR}/03_create_CNVR/cnv_clean.txt` (generated in "create CNVR" step)
  - `sample_QC.txt -> ${WKDIR}/01_initial_call/run_PennCNV/results/CNV.PennCNV_qc_new.txt` (generated when finishing PennCNV analysis; the columns "LRR_mean" and "LRR_sd" are used in this step)
  - `duplicate_pairs.txt -> ${WKDIR}/data/duplicate_pairs.txt` (optional) (tab-delimited table of two columns with header names: "sample1.name" and "sample2.name"; each row is a duplicated pair with one sample ID in the first column and the other in the second column)

Running CNV genotyping in parallel is implemented in the following four steps.

(1) Split CNVRs into different batches in each chromosome.
```sh
Rscript ${WKDIR}/04_CNV_genotype/step.1.split.cnvrs.into.batches.R \
-i ${WKDIR}/03_create_CNVR/cnvr_clean.txt \  ## generated in "create CNVR" step
-o ${WKDIR}/04_CNV_genotype/data/cnvr_batch.txt \
-n 200
```
The parameter `-n 200` indicates the maximum number of CNVRs in each batch. The script goes over the table of CNVRs in `cnvr_clean.txt` generated in the previous "create CNVR" step, appends to the table an additional column indicating the batches each CNVR belongs to, and writes the updated table to tab-delimited file `cnvr_batch.txt`.

(2) Submit parallelized jobs for CNV genotyping, each corresponding to one batch.

Note: In `step.2.submit.jobs.R`, The scripts embraced by "##<<<... ##>>>..." need to be specified based on your system.

```sh
Rscript ${WKDIR}/04_CNV_genotype/step.2.submit.jobs.R \
--type 0 \                                        ## "0" indicates initial submission
--script     ${WKDIR}/04_CNV_genotype \           ## path to main script CNV.genotype.one.chr.one.batch.R
--sourcefile ${WKDIR}/04_CNV_genotype/scripts \   ## path to relavent R functions used by the main script
--datapath   ${WKDIR}/04_CNV_genotype/data \      ## the above input files are all placed in this folder
--matrixpath ${WKDIR}/01_initial_call/finalreport_to_matrix_LRR_and_BAF/RDS \  ## LRR and BAF matrices generated in the initial step
--resultpath ${WKDIR}/04_CNV_genotype/results \   ## directory to save results
--joblog     ${WKDIR}/04_CNV_genotype/results \   ## where jobs log files to be placed
--duplicates \                                    ## (optional) indicates whether the information duplicate pairs is used in diagnosis plots
--plot                                            ## (optional) indicates whether diagnosis plots to be generated
```

When this step is finished, several subdirectories are expected to be generated:
  - `pred` (predicted copy number (CN) and associated GQ score for each individual at each CNVR)
  - `pars` (parameters of local model at each CNVR)
  - `stats` (data used to fit local model for each CNVR)
  - `png` (diagnosis plots generated in "png" format)
  - `job` (log files for submitted jobs)
  - `cnvrs_error` (list of CNVRs encourtering errors when fitting local model)

(3) Check submitted jobs and resubmit failed jobs.

Note: In `step.3.check.and.resubmit.jobs.R`, The scripts embraced by "##<<<... ##>>>..." need to be specified based on your system.

```sh
Rscript ${WKDIR}/04_CNV_genotype/step.3.check.and.resubmit.jobs.R \
--flag 1 \                                       ## 0: only print the status of submitted jobs; 1: resubmit failed jobs
--script     ${WKDIR}/04_CNV_genotype \          ## path to main script CNV.genotype.one.chr.one.batch.R
--sourcefile ${WKDIR}/04_CNV_genotype/scripts \  ## path to relavent R functions used by the main script
--datapath   ${WKDIR}/04_CNV_genotype/data \     ## the above input files are all placed in this folder
--matrixpath ${WKDIR}/01_initial_call/finalreport_to_matrix_LRR_and_BAF/RDS \  ## LRR and BAF matrices generated in the initial step
--resultpath ${WKDIR}/04_CNV_genotype/results \  ## directory to save results
--joblog     ${WKDIR}/04_CNV_genotype/results \  ## where jobs log files to be placed
--duplicates \                                   ## (optional) indicates whether the information duplicate pairs is used in diagnosis plots
--plot                                           ## (optional) indicates whether diagnosis plots to be generated
```

(4) Combine results from parallelized jobs
```sh
Rscript ${WKDIR}/04_CNV_genotype/step.4.prediction.results.R \
--datapath ${WKDIR}/04_CNV_genotype/data \     ## the above input files are all placed in this folder
--resultpath ${WKDIR}/04_CNV_genotype/results  ## directory to save results
```
When this step is finished, four files are expected to be generated in the results folder:
  - `matrix_CN.rds` (matrix of predicted copy number (CN), each row corresponds to a CNVR and each column a sample)
  - `matrix_GQ.rds` (matrix of GQ scaore, each row corresponds to a CNVR and each column a sample)
  - `cnvr_genotype.txt` (an additional column `genotype` is appended to `cnvr_batch.txt`, indicating whether CNV genotyping is successfully completed for each CNVR)
  - `sample_genotype.txt` (list of sample IDs with CN genotype).

Users can decide a threshold of GQ score afterwards. A CN genotype with GQ score below the threshold can be set as missing.

We provide an [example](https://github.com/HaoKeLab/ensembleCNV/tree/master/example/example_CNV_genotype) for this step corresponding to one example CNVR.

## 5 Boundary refinement

For a CNVR with common CNV genotype (e.g., more than 5% of CNV carriers at the CNVR), the LRR signals would be highly correlated across individuals among involved probes. We can take advantage of this structure to further refine CNVR boundaries. In other words, we are able to find a sub-block of high correlations within a local correlation matrix. When the refined boundaries are different from the initial ones, the probes falling within the range of updated boundaries will change. We need to update the local likelihood model and re-do the CNV genotyping step. If several CNVRs share the exact boundaries after boundary refinement, they will be collapsed into one. More details can be found in the [manuscript](https://doi.org/10.1093/nar/gkz068).

In current implementation, `CNVR.boundary.refinement.R` the main script that performs boundary refinement for CNVRs in one chromosome at a time. It utilizes the `Rcpp` package and implements the computationally intensive part with C++ code in `refine.cpp`. Please see this [example](https://github.com/HaoKeLab/ensembleCNV/tree/master/example/example_boundary_refinement). 

Note: 
- We split the core scripts and the platform-specific workflow scripts (see below) in order to provide flexibity for users to develop their own workflow scripts specific to their particular platform, especially when a high-performance cluster is not accessible to the users. 

- We provide workflow scripts for a cluster environment, where CNVRs within different chromosomes are processed in parallel. Relevant R and C++ scripts can be found [here](https://github.com/HaoKeLab/ensembleCNV/tree/master/05_boundary_refinement). 

Before running the script below, the following files generated in previous steps need to be copied (or linked) into the `${WKDIR}/05_boundary_refinement/data` directory, and named exactly as follows:

  - `SNP_pos.txt -> ${WKDIR}/01_initial_call/finalreport_to_matrix_LRR_and_BAF/SNP_pos.txt` (prepared in the initial step; containing chromosome and position information for each probe)
  - `cnvr_genotype.txt -> ${WKDIR}/04_CNV_genotype/results/cnvr_genotype.txt` (table of CNVR information, generated in "CNV genotyping" step)
  - `matrix_CN.rds -> ${WKDIR}/04_CNV_genotype/results/matrix_CN.rds` (matrix of CN genotype with rows as CNVRs and columns as samples, generated in "CNV genotyping" step)
  - `matrix_GQ.rds -> ${WKDIR}/04_CNV_genotype/results/matrix_GQ.rds` (matrix of GQ score with rows as CNVRs and columns as samples, generated in "CNV genotyping" step) 

Running boundary refinement in parallel is implemented in the following four steps. 

(1) Select CNVRs with common CNV genotype to be refined.
```sh
Rscript ${WKDIR}/05_boundary_refinement/step.1.common.CNVR.to.refine.R \
--datapath ${WKDIR}/05_boundary_refinement/data \       ## the above input files are all placed in this folder
--resultpath ${WKDIR}/05_boundary_refinement/results \  ## directory to save results
--freq 0.05                                             ## frequency cut-off based on which common CNVRs will be selected
```
The parameter `--freq 0.05` indicates the frequency cut-off, based on which CNVRs with common CNV genotype will be selected and subject to boundary refinement. The script goes over the table of CNVRs in `cnvr_genotype.txt` generated in the previous "CNV genotyping" step, calculates frequency of CNV genotype based on data from `matrix_CN.rds`, and appends to the table an additional column indicating the frequency of CNV genotype for each CNVR. The table of CNVRs with frequency below the cut-off will be saved in tab-delimited file `cnvr_keep.txt`, while those with frequency above the cut-off will be saved in `cnvr_refine.txt`, both in the `${WKDIR}/05_boundary_refinement/results` directory.

(2) Submit parallelized jobs for boundary refinement, each corresponding to CNVRs in one chromosome.

Note: In `step.2.submit.jobs.R`, The scripts embraced by "##<<<... ##>>>..." need to be specified based on your system.

```sh
Rscript ${WKDIR}/05_boundary_refinement/step.2.submit.jobs.R \
--refinescript ${WKDIR}/05_boundary_refinement/CNVR.boundary.refinement.R \    ## the main script for boundary refinement
--rcppfile   ${WKDIR}/05_boundary_refinement/refine.cpp \                      ## the C++ code for sub-block searching in local correlation matrix
--datapath   ${WKDIR}/05_boundary_refinement/data \                            ## the above input files are all placed in this folder
--matrixpath ${WKDIR}/01_initial_call/finalreport_to_matrix_LRR_and_BAF/RDS \  ## LRR and BAF matrices generated in the initial step
--centromere ${WKDIR}/data/centromere_hg19.txt \                               ## for other assemblies, check UCSC genome browser (see above)
--resultpath ${WKDIR}/05_boundary_refinement/results \                         ## directory to save results
--plot                                                                         ## (optional) indicates whether diagnosis plots to be generated
```

When this step is finished, several subdirectories are expected to be generated in each `res_refine/chr*` directory:
  - `data` (refined boundary information for each CNVR)
  - `png` (diagnosis plots generated in "png" format)
  - `log` (log files for submitted jobs)

(3) Combine results from parallelized jobs.
```sh
Rscript ${WKDIR}/05_boundary_refinement/step.3.clean.results.R \
--resultpath ${WKDIR}/05_boundary_refinement/results
```
When this step is finished, three files are expected to be generated in the result folder:
  - `cnvr_kept_after_refine.txt`
  - `cnvr_refined_after_refine.txt`
  - `cnvr_regenotype_after_refine.txt`

The information of rare CNVRs with CNV genotype frequency less than the cut-off will be saved in `cnvr_kept_after_refine.txt`. Some of common CNVRs go through the boundary refinement procedure, but their boundaries may remain the same, suggesting the boundaries constructed from the "create CNVR" step are correct and do not need to be updated -- the information of common CNVRs falling in this category will also be saved in `cnvr_kept_after_refine.txt`. Among CNVRs with updated boundaries, two further checking steps will be performed: a) If several CNVRs share the exact boundaries, they will be collapsed into one; b) If a CNVR has the exact boundaries as one in `cnvr_kept_after_refine.txt`, it will be removed from the list. As a result, the total number of CNVRs may be slightly reduced after boundary refinement and CNV regenotyping steps. The remaining CNVRs will be saved in `cnvr_refined_after_refine.txt`, from which a new table in the similar format as `cnvr_clean.txt` (generated in "create CNVR" step) will be generated with updated boundary information in columns `posStart`, `posEnd`, `snp_start` and `snp_end`, and be saved in `cnvr_regenotype_after_refine.txt`. The CNVRs listed in `cnvr_regenotype_after_refine.txt` will need to go through all the steps in [CNV genotyping](#4-cnv-genotyping-for-each-cnvr) to update relevant CN and GQ matrices, as the probes involved in each of these CNVRs are altered after boundary refinement. Assume the working folder for CNVR regenotyping after boundary refinement is `${WKDIR}/05a_regenotype_after_refinement`, which has similar structure as `${WKDIR}/04_CNV_genotype`. The updated CN and GQ matrices will be generated at `${WKDIR}/05a_regenotype_after_refinement/results`.

(4) Update CN and GQ matrices as well as CNVR information.
```sh
Rscript ${WKDIR}/05_boundary_refinement/step.4.update.genotype.matrix.R \
--matrixbeforerefine ${WKDIR}/05_boundary_refinement/data \        ## matrix_CN.rds and matrix_GQ.rds before boundary refinement have been copied (or linked) here
--matrixrefine ${WKDIR}/05a_regenotype_after_refinement/results \  ## path to updated CN and GQ matrices for CNVRs listed in cnvr_regenotype_after_refine.txt
--refinepath ${WKDIR}/05_boundary_refinement/results \             ## where cnvr_kept_after_refine.txt is located
--output ${WKDIR}/05_boundary_refinement/results                   ## path to save final results
```

When this step is finished, three files are expected to be generated in the results folder:
  - `matrix_CN_final.rds`
  - `matrix_GQ_final.rds`
  - `cnvr_final.txt`

In this step, the subset of CN and GQ matrices for CNVRs in `cnvr_kept_after_refine.txt` will be extracted, and combined with those generated for the CNVRs subject to re-genotyping (see above). The combined matrices are saved in `matrix_CN_final.rds` and `matrix_GQ_final.rds`. The information for the combined list of CNVRs is saved in `cnvr_final.txt`.

We provide an [example](https://github.com/HaoKeLab/ensembleCNV/tree/master/example/example_boundary_refinement) for this step at one example CNVR.

## 6 Performance assessment

In large-scale genetic studies, for QC purposes in SNP genotyping, technical duplicates are often available. The concordance rate of CNV calls in duplicated pairs, i.e. the reproducibility, can be used as a surrogate of accuracy measurement. In ensembleCNV, we define the genotyping quality (GQ) score to quantify the confidence of the CN genotype assigned to each individual at each CNVR. As the GQ score threshold increases, the concordance rate constantly increases at the cost of decreased sample-wise and CNVR-wise call rates. The users can select a GQ score threshold to achieve a balance between concordance rate and call rate. More details can be found in the [manuscript](https://doi.org/10.1093/nar/gkz068).

When technical duplicates are available, the users can use the following script to evaluate the quality of CNV calls.

(1) Evaluate concordance rate of CNV calls between technical duplicates as well as sample-wise and CNVR-wise call rates.
```sh
Rscript ${WKDIR}/06_performance_assessment/step.1.performance.assessment.R \
--duplicates ${WKDIR}/data/duplicate_pairs.txt \                          ## duplicates information, an option in "CNV genotyping" step (see above)
--matrixCN ${WKDIR}/05_boundary_refinement/results/matrix_CN_final.rds \  ## CN matrix generated in "boundary refinement" step (see above)
--matrixGQ ${WKDIR}/05_boundary_refinement/results/matrix_GQ_final.rds \  ## GQ matrix generated in "boundary refinement" step (see above)
--resultpath ${WKDIR}/06_performance_assessment                           ## path to directory for saving results
```
When this step is finished, two files will be generated in results folder:
 - `performance_assessment.rds` (concordance rate, number of CNVRs, sample-wise call rate and CNVR-wise call rate at different GQ score thresholds)
 - `performance_assessment.png` (visualization of information in `performance_assessment.rds`, based on which the users can choose a GQ score threhold to achieve desirable between concordance rate and call rate)

When technical duplicates are not available, the users can skip step (1) and choose an empirical GQ score threshold directly. 

(2) Set GQ score threshold to generate final results.
```sh
Rscript ${WKDIR}/06_performance_assessment/step.2.set.GQ.generate.results.R \
--matrixCN ${WKDIR}/05_boundary_refinement/results/matrix_CN_final.rds \  ## CN matrix generated in "boundary refinement" step
--matrixGQ ${WKDIR}/05_boundary_refinement/results/matrix_GQ_final.rds \  ## GQ matrix generated in "boundary refinement" step
--cnvrfile ${WKDIR}/05_boundary_refinement/results/cnvr_final.txt \       ## CNVR information generated in "boundary refinement" step
--resultpath ${WKDIR}/06_performance_assessment \                         ## path to directory for saving results
--gqscore <INT>  ## GQ score threhold chosen based on evaluation in step (1) or chosen empirically based on previous studies
```
When this step is finished, three files will be generated in results folder:
 - `matrix_CN_after_GQ.rds` (CN matrix with rows as CNVRs and columns as samples; the CN values associated with GQ < threshold are set as missing value, which is indicated by -9)
 - `cnvr_after_GQ.txt` (CNVR information and summary statistics; CNVRs with no CNV calls (i.e. no CN = 0, 1, 3) are excluded)
 - `sample_after_GQ.txt`(sample information and summary statistics)

## Example

We provide examples for three major parts of the pipeline: [CNVR creating](https://github.com/HaoKeLab/ensembleCNV/tree/master/example/example_create_CNVR), [CNV genotyping](https://github.com/HaoKeLab/ensembleCNV/tree/master/example/example_CNV_genotype), and [boundary refinement](https://github.com/HaoKeLab/ensembleCNV/tree/master/example/example_boundary_refinement). 

