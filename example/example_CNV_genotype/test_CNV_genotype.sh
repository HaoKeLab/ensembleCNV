## path to the directory in which all ensembleCNV scripts and data are located. 
path_ensembleCNV=""

Rscript ${path_ensembleCNV}/04_CNV_genotype/CNV.genotype.one.chr.one.batch.R \
--chr 1 \
--batch 2 \
--type 0 \
--datapath ${path_ensembleCNV}/example/example_CNV_genotype/data/ \
--resultpath ${path_ensembleCNV}/example/example_CNV_genotype/results/ \
--matrixpath ${path_ensembleCNV}/example/example_CNV_genotype/matrix/ \
--sourcefile ${path_ensembleCNV}/04_CNV_genotype/scripts/ \
--duplicates \
--plot
