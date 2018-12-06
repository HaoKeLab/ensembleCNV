## Example: Boundary refinement

Here is a demo of the main script `CNVR.boundary.refinement.R` for [boundary refinement](https://github.com/HaoKeLab/ensembleCNV#5-boundary-refinement) using example data from one CNVR.

Please specify where the git clone of ensembleCNV is located.
```sh
ENSEMBLECNV=</path/to/ensembleCNV>
```

Then run the following code for a demo of the main script `CNVR.boundary.refinement.R` for boundary refinement.
```sh
Rscript ${ENSEMBLECNV}/05_boundary_refinement/CNVR.boundary.refinement.R \
--chr 2 \
--rcppfile ${ENSEMBLECNV}/05_boundary_refinement/refine.cpp \
--datapath ${ENSEMBLECNV}/example/example_boundary_refinement/data \
--matrixpath ${ENSEMBLECNV}/example/example_boundary_refinement/RDS \
--centromere ${ENSEMBLECNV}/example/example_boundary_refinement/data/centromere_hg19.txt \
--resultpath ${ENSEMBLECNV}/example/example_boundary_refinement/results \
--plot
```

Note: 

- When the analysis is successfully completed, the output will be stored at the directory `${ENSEMBLECNV}/example/example_boundary_refinement/results/res_refine`.

- In practice, the list of common CNVRs in `cnvr_refine.txt`, whose boundaries are to be refined, is selected by the step `${ENSEMBLECNV}/05_boundary_refinement/step.1.common.CNVR.to.refine.R` based on frequency cut-off specified by the user, before boundary refinement is actually performed (see step (1) of [boundary refinement](https://github.com/HaoKeLab/ensembleCNV#5-boundary-refinement)). Therefore, `cnvr_refine.txt` is supposed to appear in the directory `${ENSEMBLECNV}/example/example_boundary_refinement/results` (instead of the `data` folder) as input for subsequent `CNVR.boundary.refinement.R`. 
