
## Example: Boundary refinement

Here is a demo of [creating CNVR](https://github.com/HaoKeLab/ensembleCNV#3-create-cnvr) using example data of CNVs clumping around one CNVR.

Please specify where the git clone of ensembleCNV is located.
```sh
ENSEMBLECNV=</path/to/ensembleCNV>
```

Step 1: reformat the CNV calls generated from individual CNV caller: iPattern, PennCNV and QuantiSNP
```sh
Rscript ${ENSEMBLECNV}/03_create_CNVR/step.1.CNV.data.R \
${ENSEMBLECNV}/example/example_create_CNVR/results/ \
${ENSEMBLECNV}/example/example_create_CNVR/data/iPattern_all_calls.txt \
${ENSEMBLECNV}/example/example_create_CNVR/data/CNV.PennCNV_new.txt \
${ENSEMBLECNV}/example/example_create_CNVR/data/quantisnp.cnv \
${ENSEMBLECNV}/example/example_create_CNVR/data/Samples_Table.txt
```
Note:

- `iPattern_all_calls.txt`, `CNV.PennCNV_new.txt`, and `quantisnp.cnv` are examples of what raw CNV calls generated from iPattern, PennCNV and QuantiSNP look like.

- We do not include `Gender` column in `Samples_Table.txt` as gender information is not relevant for creating CNVR in this example.

- When this step is successfully completed, you will find in `${ENSEMBLECNV}/example/example_create_CNVR/results/` directory `cnv.ipattern.txt`, `cnv.penncnv.txt`, and `cnv.quantisnp.txt`, which are reformated CNV calls from the 3 CNV callers.

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

Note: When this step is successfully completed, you will find some intermediate outputs and the final results in the directory `${ENSEMBLECNV}/example/example_create_CNVR/results/`, including
- `cnv_clean.txt`: the table of merged CNV events from iPattern, PennCNV and QuantiSNP; the `CNVR_ID` in the table indicates which CNVR each CNV belongs to.
- `cnvr_clean.txt`: the table of constructed CNVRs with each assigned a `CNVR_ID`.
