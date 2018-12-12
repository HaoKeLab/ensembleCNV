## QuantiSNP

### Installation

To download and install QuantiSNP, please follow the detailed instructions at the [page](https://sites.google.com/site/quantisnp/downloads), which provides links to MATLAB Run-Time Component (MCR) Libraries, QuantiSNP package and GC content data. For more information about QuantiSNP, please refer to their original [QuantiSNP website](https://sites.google.com/site/quantisnp/home).

After installation, set up environment variable QUANTISNP: `export QUANTISNP='/path/to/quantisnp'`

Please organize the installation folder in the following way:

- MATLAB Run-Time Component (MCR) Libraries root directory: `${QUANTISNP}/v79/`
- QuantiSNP root directory: `${QUANTISNP}/quantisnp/`
- GC content data (take b37/hg19 data for example) directory: `${QUANTISNP}/data/b37`

Note:

- Running QuantiSNP does not require MATLAB, but rather the developers provided a self-contained MATLAB Run-Time Component Libraries in accompany with QuantiSNP.

- We have successfully installed MCR and QuantiSNP on the Linux CentOS 6.9 with openjdk 1.6 or Ubuntu 16.04 with openjdk 1.8.

### Analysis workflow

Note: The auxiliary scripts we provide here were used on our high performance cluster. The users need to modifiy the scripts according to the specific system the users are using. Please refer to original [QuantiSNP website](https://sites.google.com/site/quantisnp/) for more information.

Note: 

- PennCNV was originally designed to sequentially analyze one sample at a time. Please refer to [PennCNV website](http://penncnv.openbioinformatics.org/en/latest/) for how to run PennCNV in a sequential way. Here, we provide scripts to run the analysis on multiple samples in parallel via job submitting system (one sample per job) in a cluster environment. 

- In the following steps (1) and (2), the scripts regarding job submission embraced by "##<<<... ##>>>..." in the scripts need to be specified by the users based on the system the users are using.

We run PennCNV analysis with the following 5 steps:


Running QuantiSNP includes the following 3 steps:

(1) Run QuantiSNP for each sample in parallel (through job scheduling system on cluster)
```sh
Rscript ${WKDIR}/01_initial_call/run_QuantiSNP/step.1.prepare.QuantiSNP.R \
--quantisnp ${QUANTISNP} \
--data ${WKDIR}/01_initial_call/run_QuantiSNP/data \ ## generated with finalreport_to_QuantiSNP.pl
--sample ${WKDIR}/data/Samples_Table.txt \
--result ${WKDIR}/01_initial_call/run_QuantiSNP/results/res
```
(2) Check job status and resubmit unfinishing jobs
```sh
Rscript ${WKDIR}/01_initial_call/run_QuantiSNP/step.2.check.QuantiSNP.R \
--quantisnp ${QUANTISNP} \
--data ${WKDIR}/01_initial_call/run_QuantiSNP/data \ ## generated with finalreport_to_QuantiSNP.pl
--sample ${WKDIR}/data/Samples_Table.txt \
--result ${WKDIR}/01_initial_call/run_QuantiSNP/results/res
```

(3) Combine PennCNV results from each sample, including the content in ".cnv" files
```sh
perl ${WKDIR}/01_initial_call/run_QuantiSNP/step.3.combine.QuantiSNP.pl \
--in_dir ${WKDIR}/01_initial_call/run_QuantiSNP/results/res \
--out_dir ${WKDIR}/01_initial_call/run_QuantiSNP/results
```
