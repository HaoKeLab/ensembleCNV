### Run QuantiSNP

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
