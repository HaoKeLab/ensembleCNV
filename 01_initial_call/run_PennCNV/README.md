### Run PennCNV

Note: The auxiliary scripts we provide here were used on our high performance cluster. The users need to modifiy the scripts according to the specific system the users are using. Please refer to original [PennCNV website](http://penncnv.openbioinformatics.org/en/latest/) for more information.

Running PennCNV includes the following 5 steps:

(1) Prepare SNP.pfb and SNP.gcmodel files

See scripts in `step.1.prepare.files.sh` for details.

(2) Run PennCNV for each sample in parallel (through job scheduling system on cluster)
```sh 
Rscript step.2.run.PennCNV.jobs.R \
--data /path_to_data/ \ ## generated with finalreport_to_PennCNV.pl
--wkdir /wk_dir/ \
--pfb /path/to/SNP.pfb \
--gcmodel /path/to/SNP.gcmodel \
--hmm /path_to_penncnv/lib/hhall.hmm
```

(3) Check job status and resubmit unfinishing jobs
```sh
Rscrip step.3.check.PennCNV.jobs.R \
--data /path_to_data/ \ ## generated with finalreport_to_PennCNV.pl
--wkdir /wk_dir/ \
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
