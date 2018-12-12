## QuantiSNP

### Installation

To download and install QuantiSNP (version 2), please follow the detailed instructions at the [page](https://sites.google.com/site/quantisnp/downloads), which provides links to download MATLAB Run-Time Component Libraries, QuantiSNP package and GC content data. For more information about QuantiSNP, please refer to their original [QuantiSNP website](https://sites.google.com/site/quantisnp/home).

After installation, set up environment variable QUANTISNP: `export QUANTISNP='/path/to/quantisnp'`

Please organize the installation folder in the following way:

- MATLAB Run-Time Component Libraries root directory: `${QUANTISNP}/v79/`
- QuantiSNP root directory: `${QUANTISNP}/quantisnp/`
- GC content data (take b37/hg19 data for example) directory: `${QUANTISNP}/data/b37/`

Note:

- Running QuantiSNP does not require MATLAB, but rather the developers provided a self-contained MATLAB Run-Time Component Libraries in accompany with QuantiSNP.

- We have checked that the installation of MATLAB Run-Time Component Libraries and QuantiSNP worked properly on two versions of Linux: CentOS 6.9 with openjdk 1.6 (the system used on [Minverva](https://hpc.mssm.edu/) cluster) or Ubuntu 16.04 with openjdk 1.8. The installation of the two components will probably require some further tweaking for different systems, due to tricky JRE (Java Runtime Environment) and libxp dependencies.

### Analysis workflow

Note: 

- QuantiSNP was originally designed to analyze one sample at a time or a batch of samples sequentially. Please refer to the original QuantiSNP [usage](https://sites.google.com/site/quantisnp/howto) for more details. Here, we provide scripts to run the analysis on multiple samples in parallel via job submitting system (one sample per job) in a cluster environment. 

- In the following steps (1) and (2), the scripts regarding job submission embraced by "##<<<... ##>>>..." in the scripts need to be specified by the users based on the system the users are using.

We run QuantiSNP analysis with the following 3 steps:

(1) Run QuantiSNP for each sample in parallel (through job scheduling system on cluster)
```sh
Rscript ${WKDIR}/01_initial_call/run_QuantiSNP/step.1.prepare.QuantiSNP.R \
--quantisnp ${QUANTISNP} \
--data ${WKDIR}/01_initial_call/run_QuantiSNP/data \ ## generated with finalreport_to_QuantiSNP.pl
--sample ${WKDIR}/data/Samples_Table.txt \
--result ${WKDIR}/01_initial_call/run_QuantiSNP/results/res
```
When the analysis is completed, there will be subfolders named after sample IDs, each for one sample respectively, created in the directory `${WKDIR}/01_initial_call/run_QuantiSNP/results/res`. Within each sample subfolders, two files (among others) will be generated and used in downstream analysis:
- `<Sample_ID>.qc`: chromosome-level summary statistics, which will be summarized later at sample level and used in checking [batch effect](https://github.com/HaoKeLab/ensembleCNV#pca-on-summary-statistics). 
- `<Sample_ID>.cnv`: raw CNV calls for each sample.

(2) Check job status and resubmit unfinishing jobs
```sh
Rscript ${WKDIR}/01_initial_call/run_QuantiSNP/step.2.check.QuantiSNP.R \
--quantisnp ${QUANTISNP} \
--data ${WKDIR}/01_initial_call/run_QuantiSNP/data \ ## generated with finalreport_to_QuantiSNP.pl
--sample ${WKDIR}/data/Samples_Table.txt \
--result ${WKDIR}/01_initial_call/run_QuantiSNP/results/res
```
This step checks if the jobs submitted for each sample in step (1) are successfully completed and resubmits failed jobs if there is any.

(3) Combine PennCNV results from each sample, including the content in ".cnv" files
```sh
perl ${WKDIR}/01_initial_call/run_QuantiSNP/step.3.combine.QuantiSNP.pl \
--in_dir ${WKDIR}/01_initial_call/run_QuantiSNP/results/res \
--out_dir ${WKDIR}/01_initial_call/run_QuantiSNP/results
```
When the analysis is completed, you will find `quantisnp.cnv`, which will be used by ensembleCNV, in the directory `${WKDIR}/01_initial_call/run_QuantiSNP/results`. `quantisnp.cnv` combines the CNV calls from all samples generated in steps (1) and (2).
