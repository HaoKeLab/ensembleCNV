## path to the directory in which all ensembleCNV scripts and data are located. 
path_ensembleCNV=""

Rscript ${path_ensembleCNV}/05_boundary_refinement/CNVR.boundary.refinement.R \
--chr 2 \
--datapath ${path_ensembleCNV}/example/example_boundary_refinement/data/ \
--matrixpath ${path_ensembleCNV}/example/example_boundary_refinement/matrix/ \
--resultpath ${path_ensembleCNV}/example/example_boundary_refinement/results/ \
--rcppfile ${path_ensembleCNV}/05_boundary_refinement/refine.cpp \
--centromere ${path_ensembleCNV}/example/example_boundary_refinement/data/centromere_hg19.txt \
--plot
