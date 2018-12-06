
## Example: Boundary refinement

Here is a demo of [creating CNVR]() using example data of CNVs clumping around one CNVR.

Please specify where the git clone of ensembleCNV are located.
```sh
ENSEMBLECNV=</path/to/ensembleCNV>
```


Rscript ${path_ensembleCNV}/03_create_CNVR/step.1.CNV.data.R \
${path_ensembleCNV}/example/example_create_CNVR/results/ \
${path_ensembleCNV}/example/example_create_CNVR/data/iPattern_cnvr1_all_calls.txt \
${path_ensembleCNV}/example/example_create_CNVR/data/PennCNV_cnvr1.txt \
${path_ensembleCNV}/example/example_create_CNVR/data/QuantiSNP_cnvr1.txt \
${path_ensembleCNV}/example/example_create_CNVR/data/Sample_Map.txt

Rscript ${path_ensembleCNV}/03_create_CNVR/step.2.create.CNVR.R \
--icnv ${path_ensembleCNV}/example/example_create_CNVR/results/cnv.ipattern.txt \
--pcnv ${path_ensembleCNV}/example/example_create_CNVR/results/cnv.penncnv.txt \
--qcnv ${path_ensembleCNV}/example/example_create_CNVR/results/cnv.quantisnp.txt \
--snp ${path_ensembleCNV}/example/example_create_CNVR/data/SNP.chr1.pfb \
--centromere ${path_ensembleCNV}/example/example_create_CNVR/data/centromere_hg19.txt \
--output ${path_ensembleCNV}/example/example_create_CNVR/results/