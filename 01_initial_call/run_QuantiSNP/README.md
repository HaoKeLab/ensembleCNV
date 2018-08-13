### Run QuantiSNP

Note: The auxiliary scripts we provide here were used on our high performance cluster. The users need to modifiy the scripts according to the specific system the users are using. Please refer to original [QuantiSNP website](https://sites.google.com/site/quantisnp/) for more information.

Running QuantiSNP includes the following 3 steps:

(1) Run PennCNV for each sample in parallel (through job scheduling system on cluster)
```sh
Rscript step.1.prepare.QuantiSNP.R \
-i /path/to/data/ \ ## generated with finalreport_to_QuantiSNP.pl
-o /path/to/results/
```
(2) Check job status and resubmit unfinishing jobs
```sh
Rscript step.2.check.QuantiSNP.R \
-d /path/to/data/ \ ## generated with finalreport_to_QuantiSNP.pl
-r /path/to/results/ 
```

(3) Combine PennCNV results from each sample, including the content in .cnv files
```sh
## To run this script, you need to add "in_dir" and "out_dir" information in the script.
perl step.3.combine.QuantiSNP.pl
```
