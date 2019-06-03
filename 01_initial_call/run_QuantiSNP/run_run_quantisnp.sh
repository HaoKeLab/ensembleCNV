#!/bin/bash 

$QUANTISNP/quantisnp/linux64/run_quantisnp2.sh $MCR --chr $CHR --outdir $OUT --sampleid $SAMPLEID --gender $GENDER --emiters 10 --lssetting 2000000 --gcdir $GC --plot --genotype --config $PARAM --levels $LEV --input-files $IN --chrX 23 --doXcorrect
