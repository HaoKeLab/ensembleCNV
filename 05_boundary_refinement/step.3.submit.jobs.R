

script_refinement <- ""

## input for script_refinement
file_cnvrs_refine <- ""
path_lrr <- ""
file_pfb <- ""
file_centromere = ""
path_png <- ""
path_result <- ""
file_rcpp <- ""

for ( chr1 in 1:22 ) {
  
  cmd.chr1 <- paste( script_refinement, 
                     "-c", chr1,
                     "-r", file_cnvrs_refine,
                     "-l", path_lrr,
                     "-p", file_pfb,
                     "-m", file_centromere,
                     "-g", path_png,
                     "-o", path_result,
                     "-s", file_rcpp)
  
  bsub.cmd <- paste("bsub -n 2 -W 10:00 -R 'rusage[mem=10000]' -P acc_haok01a", 
                    "-J", chr1,
                    "-q premium",
                    shQuote( cmd.chr1 ))
  
  cat( bsub.cmd, "\n" )
  system( bsub.cmd )
  Sys.sleep(0.5)
} 






