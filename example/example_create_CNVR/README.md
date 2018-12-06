
## Example: Boundary refinement

Here is a demo of [creating CNVR](https://github.com/HaoKeLab/ensembleCNV#3-create-cnvr) using example data of CNVs clumping around one CNVR.

Please specify where the git clone of ensembleCNV are located.
```sh
ENSEMBLECNV=</path/to/ensembleCNV>
```

Step 1: reformat the CNV calls generated from individual CNV caller: iPattern, PennCNV and QuantiSNP
```sh
Rscript ${ENSEMBLECNV}/03_create_CNVR/step.1.CNV.data.R \
${ENSEMBLECNV}/example/example_create_CNVR/results/ \
${ENSEMBLECNV}/example/example_create_CNVR/data/iPattern_cnvr1_all_calls.txt \
${ENSEMBLECNV}/example/example_create_CNVR/data/PennCNV_cnvr1.txt \
${ENSEMBLECNV}/example/example_create_CNVR/data/QuantiSNP_cnvr1.txt \
${ENSEMBLECNV}/example/example_create_CNVR/data/Samples_Table.txt
```

Step 2: create CNVR
```sh
Rscript ${ENSEMBLECNV}/03_create_CNVR/step.2.create.CNVR.R \
--icnv ${ENSEMBLECNV}/example/example_create_CNVR/results/cnv.ipattern.txt \
--pcnv ${ENSEMBLECNV}/example/example_create_CNVR/results/cnv.penncnv.txt \
--qcnv ${ENSEMBLECNV}/example/example_create_CNVR/results/cnv.quantisnp.txt \
--snp ${ENSEMBLECNV}/example/example_create_CNVR/data/SNP_pos.txt \
--centromere ${ENSEMBLECNV}/example/example_create_CNVR/data/centromere_hg19.txt \
--output ${ENSEMBLECNV}/example/example_create_CNVR/results/
```

