#!/bin/bash

Rscript CNV.genotype.one.chr.one.batch.R \
            --type $TYPE \ 
            --datapath $DTPATH \ 
            --resultpath $RESPATH \ 
            --matrixpath $MATPATH \ 
            --sourcefile $SOURCE
