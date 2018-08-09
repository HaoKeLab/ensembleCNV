### call PennCNV

Here, calling PennCNV including following 5 steps:

prepare files containing SNP.pfb and SNp.gcmodel for running PennCNV:
```sh
step.0.prepare.files.sh contains all commands 
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

combine all PennCNV calling results (sample based):
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
