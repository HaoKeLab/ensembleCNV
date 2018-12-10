## PennCNV

### Installation

To install PennCNV, please follow the detailed instructions (including trouble shooting) at the PennCNV [page](http://penncnv.openbioinformatics.org/en/latest/user-guide/install/). For more information about PennCNV, please refer to their original [PennCNV website](http://penncnv.openbioinformatics.org/en/latest/).

After installation, set up environment variable PENNCNV: `export PENNCNV='/path/to/penncnv'`

### Analysis workflow

Note: The auxiliary scripts we provide here were used on our high performance cluster. The users need to modifiy the scripts according to the specific system the users are using. 

Running PennCNV includes the following 5 steps:

(1) Prepare SNP.pfb and SNP.gcmodel files

See scripts in `step.1.prepare.files.sh` for details.


Note: PennCNV was originally designed to sequentially analyze one sample at a time. Please refer to [PennCNV website](http://penncnv.openbioinformatics.org/en/latest/) for how to perform a sequential analysis. Here, we provide scripts to run the analysis on multiple samples in parallel via job submitting system (one sample per job) in a cluster environment. 

In the following steps (2) and (3), the scripts embraced by "##<<<... ##>>>..." in the scripts need to be specified based on your system.

(2) Run PennCNV for each sample in parallel (through job submitting system on cluster)

Note: In `step.2.run.PennCNV.jobs.R`, The scripts embraced by "##<<<... ##>>>..." need to be specified based on your system.

```sh 
Rscript ${WKDIR}/01_initial_call/run_PennCNV/step.2.run.PennCNV.jobs.R \
--penncnv ${PENNCNV} \
--data ${WKDIR}/01_initial_call/run_PennCNV/data/ \     ## generated with finalreport_to_PennCNV.pl
--wkdir ${WKDIR}/01_initial_call/run_PennCNV/results/ \ ## output directory
--pfb ${WKDIR}/01_initial_call/run_PennCNV/data_aux/SNP.pfb \
--gcmodel ${WKDIR}/01_initial_call/run_PennCNV/data_aux/SNP.gcmodel \
--hmm ${PENNCNV}/lib/hhall.hmm
```

(3) Check job status and resubmit failed jobs

Note: In `step.3.check.PennCNV.jobs.R`, The scripts embraced by "##<<<... ##>>>..." need to be specified based on your system.

```sh
Rscrip step.3.check.PennCNV.jobs.R \
--data /path_to_data/ \ ## generated with finalreport_to_PennCNV.pl
--wkdir /wk_dir/ \ ## output directory
--pfb /path/to/SNP.pfb \
--gcmodel /path/to/SNP.gcmodel \
--hmm /path_to_penncnv/lib/hhall.hmm
```

(4) Combine PennCNV results from each sample, including the content in .rawcnv and .log files
```sh
perl step.4.combine.PennCNV.res.pl \
--in_dir /path/to/results/ \
--out_dir /path/to/output/
```

(5) Merge closely adjacent CNVs and generate final results
```sh
Rscript step.5.clean.PennCNV.res.R \
-i /path/to/results/ \
-f /path/to/SNP.pfb \
-n file_name_of_combined_results ## [file_name_of_combined_results].rawcnv from step 4
```
