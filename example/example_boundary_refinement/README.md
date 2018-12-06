## Example: Boundary refinement

Please specify where the scripts of ensembleCNV are located.
```sh
ENSEMBLECNV=</path/to/ensembleCNV>
```

This is a demo of how the main script `CNVR.boundary.refinement.R` for boundary refinement works.
```sh
Rscript ${ENSEMBLECNV}/05_boundary_refinement/CNVR.boundary.refinement.R \
--chr 2 \
--rcppfile ${ENSEMBLECNV}/05_boundary_refinement/refine.cpp \
--datapath ${ENSEMBLECNV}/example/example_boundary_refinement/data/ \
--matrixpath ${ENSEMBLECNV}/example/example_boundary_refinement/matrix/ \
--centromere ${ENSEMBLECNV}/example/example_boundary_refinement/data/centromere_hg19.txt \
--resultpath ${ENSEMBLECNV}/example/example_boundary_refinement/results/ \
--plot
```
