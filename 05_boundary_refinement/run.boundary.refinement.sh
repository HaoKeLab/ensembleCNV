#!/bin/bash

Rscript ${WKDIR}/05_boundary_refinement/CNVR.boundary.refinement.R \ 
            --datapath $DTPATH \ 
            --resultpath $RESPATH \ 
            --matrixpath $MATPATH \ 
            --rcppfile $RCPPFILE \ 
            --centromere $CENTROMERE 
