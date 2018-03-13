#!/hpc/packages/minerva-common/R/3.3.1/lib64/R/bin/Rscript --vanilla

path_code   ="/sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/ensembleCNV/pipeline_regenotypeCNV_raw"
name_script = "runme_minerva_each_chr_each_batch.R"

script = file.path(path_code, name_script)

path_cnvr = "/sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/ensembleCNV/5batch/res"
file_cnvr = "cnvrs_batch_annotate.rds"  # with batch annotated

dat_cnvr = readRDS(file = file.path(path_cnvr, file_cnvr)) 

chrs = sort( unique(dat_cnvr$chr) )

for ( chr1 in chrs ) {
  
  dat_cnvr_chr1 = subset(dat_cnvr, chr == chr1)
  batch_chr1 = sort( unique(dat_cnvr_chr1$batch) )
  
  if ( nrow(dat_cnvr_chr1) == 0) {
    next
  }
  
  for ( batch1 in batch_chr1 ){
    
    cat("chr:", chr1, "batch1:", batch1, "\n")
    cmd1 = paste(script, "-c", chr1, "-b", batch1, "-t", 0)
    bsub.cmd = paste("bsub -n 2 -W 10:00 -R 'rusage[mem=20000]' -P acc_haok01a",
                     "-q alloc", shQuote(cmd1))
    cat(bsub.cmd, "\n")
    system(bsub.cmd)
    
    Sys.sleep(0.1)
  }
  
}

# for test -----------------------------------------------------
# path_code <- "/sc/orga/projects/haok01a/chengh04/Food_Allergy.New/Code_CNV_gatk"
# file_script <- "runme_minerva_each_chr_each_batch.R"
# 
# script <- paste(path_code, file_script, sep = "/")
# cmd <- paste(script, 1, 1)
# bsub.cmd <- paste("bsub -n 2 -W 20:00 -R 'rusage[mem=10000]' -P acc_haok01a",
#                   "-q alloc", shQuote(cmd))
# system(bsub.cmd)



