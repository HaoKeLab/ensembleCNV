## Example: CNV genotyping

Here is a demo of the main script `CNV.genotype.one.chr.one.batch.R` for [CNV genotyping](https://github.com/HaoKeLab/ensembleCNV#4-cnv-genotyping-for-each-cnvr) using example data from one CNVR.

Please specify where the git clone of ensembleCNV is located.
```sh
ENSEMBLECNV=</path/to/ensembleCNV>
```

Then run the following code for a demo of the main script `CNV.genotype.one.chr.one.batch.R` for CNV genotyping.
```sh
Rscript ${path_ensembleCNV}/04_CNV_genotype/CNV.genotype.one.chr.one.batch.R \
--chr 1 \
--batch 2 \
--type 0 \
--sourcefile ${path_ensembleCNV}/04_CNV_genotype/scripts/ \
--datapath ${path_ensembleCNV}/example/example_CNV_genotype/data \
--matrixpath ${path_ensembleCNV}/example/example_CNV_genotype/RDS \
--resultpath ${path_ensembleCNV}/example/example_CNV_genotype/results \
--duplicates \
--plot
```

Note: When the analysis is successfully completed, in the directory `${path_ensembleCNV}/example/example_CNV_genotype/results`, you will find similar directory structure as in a real project.
