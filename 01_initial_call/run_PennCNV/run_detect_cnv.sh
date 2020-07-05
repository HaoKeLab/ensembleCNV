#!/bin/bash

perl $PENNCNV/detect_cnv.pl -test --confidence -hmm $HMM -pfb $PFB -gcmodel $GC -list $LIST -log $LOG -out $OUT
