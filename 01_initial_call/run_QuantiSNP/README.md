### Run QuantiSNP

Note: The auxiliary scripts we provide here were used on our high performance cluster. The users need to modifiy the scripts according to the specific system the users are using. Please refer to original [QuantiSNP website](https://sites.google.com/site/quantisnp/) for more information.

Running QuantiSNP includes the following 3 steps:

(1) Prepare QuantiSNP and submit jobs:
```sh
Rscript step.1.prepare.QuantiSNP.R \
-i /path/to/data/ \ ## generated with finalreport_to_QuantiSNP.pl
-o /path/to/results/
```
(2) Check jobs and resubmit:
```sh
Rscript step.2.check.QuantiSNP.R \
-d /path/to/data/ \
-r /path/to/callingCNV/folder 
```

(3) Combine CNV calling results:
running this script, you need to add "in_dir", "out_dir", "out_file" information in the script.
```sh
perl step.3.combine.QuantiSNP.pl
```
