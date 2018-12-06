## Example: Boundary refinement

Please specify where the scripts of ensembleCNV are located.
```sh
ENSEMBLECNV=</path/to/ensembleCNV>
```

Then run the following code for a demo of the main script `CNVR.boundary.refinement.R` for boundary refinement.
```sh
Rscript ${ENSEMBLECNV}/05_boundary_refinement/CNVR.boundary.refinement.R \
--chr 2 \
--rcppfile ${ENSEMBLECNV}/05_boundary_refinement/refine.cpp \
--datapath ${ENSEMBLECNV}/example/example_boundary_refinement/data \
--matrixpath ${ENSEMBLECNV}/example/example_boundary_refinement/matrix \
--centromere ${ENSEMBLECNV}/example/example_boundary_refinement/data/centromere_hg19.txt \
--resultpath ${ENSEMBLECNV}/example/example_boundary_refinement/results \
--plot
```

Note: 

- In practice, the list of common CNVRs in `cnvr_refine.txt`, whose boundaries are to be refined, is selected by the step `${ENSEMBLECNV}/05_boundary_refinement/step.1.common.CNVR.to.refine.R` with frequency cut-off (e.g. 0.05) before boundary refinement is actually performed (see step (1) of [boundary refinement](https://github.com/HaoKeLab/ensembleCNV#5-boundary-refinement)). Therefore, `cnvr_refine.txt` is supposed to appear in the directory `${ENSEMBLECNV}/example/example_boundary_refinement/results` (instead of the `data` folder) as the input for `CNVR.boundary.refinement.R`. 

- The actual output of the above code is stored at the directory `${ENSEMBLECNV}/example/example_boundary_refinement/results/res_refine`.
