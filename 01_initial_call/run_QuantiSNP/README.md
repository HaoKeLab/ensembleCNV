### Run QuantiSNP

Note: The auxiliary scripts we provide here were used on our high performance cluster. The users need to modifiy the scripts according to the specific system the users are using. Please refer to original [QuantiSNP documents](https://sites.google.com/site/quantisnp/) for more information.

Running QuantiSNP includes 3 steps:

(1) Prepare QuantiSNP and submit jobs:
```sh
./step.1.prepare.QuantiSNP.R \
-i path/to/data/folder \
-o path/to/result/folder
```
(2) Check jobs and resubmit:
```sh
./step.2.check_QuantiSNP.R \
-d path/to/data/folder \
-r path/to/callingCNV/folder 
```

(3) Combine CNV calling results:
running this script, you need to add "in_dir", "out_dir", "out_file" information in the script.
```sh
perl step.3.combine.QuantiSNP.pl
```
