## Example: CNV genotyping

Here is a demo of the main script `CNV.genotype.one.chr.one.batch.R` for [CNV genotyping](https://github.com/HaoKeLab/ensembleCNV#4-cnv-genotyping-for-each-cnvr) using example data from one CNVR.

Please specify where the git clone of ensembleCNV is located.
```sh
ENSEMBLECNV=</path/to/ensembleCNV>
```

Then run the following code for a demo of the main script `CNV.genotype.one.chr.one.batch.R` for CNV genotyping.
```sh
Rscript ${ENSEMBLECNV}/04_CNV_genotype/CNV.genotype.one.chr.one.batch.R \
--chr 1 \
--batch 2 \
--type 0 \
--sourcefile ${ENSEMBLECNV}/04_CNV_genotype/scripts/ \
--datapath ${ENSEMBLECNV}/example/example_CNV_genotype/data \
--matrixpath ${ENSEMBLECNV}/example/example_CNV_genotype/RDS \
--resultpath ${ENSEMBLECNV}/example/example_CNV_genotype/results \
--duplicates \
--plot
```

Note: When the analysis is successfully completed, in the directory `${path_ensembleCNV}/example/example_CNV_genotype/results`, you will find similar directory structure and outputs as in a real project. In particular,

- in the subfolders of the `pred` folder, you will find `*_pred.rds`, each corresponding to the CN genotype and GQ score for a CNVR. They are stored in `.rds` format in order to save space and improve I/O time.

- in the subfolders of the `png` folder, you will find different diagnosis plots for each CNVR.
