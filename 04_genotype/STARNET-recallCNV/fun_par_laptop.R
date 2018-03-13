

path_intput <- "C:/mssm_projects/STARNET/work_3/Predict_All_CNVR/"

# PennCNV sample LRR
samples_LRR <- read.table(file = file.path(path_intput, 'res_one_CNVR/StarNet_PennCNV_qc.txt'),
                          sep = "\t", header = TRUE, as.is = TRUE, check.names = FALSE)
samples_LRR$Sample_ID <- gsub(pattern = ".txt", replacement = "", samples_LRR$File)

# sample1.name sample2.name
dup_pairs <- read.table(file = file.path(path_intput, 'res_one_CNVR/samples_dup_pair_qc.txt'),
                        sep = "\t", header = TRUE, as.is = TRUE, check.names = FALSE)

sd_par <- 0.8
paras_LRR <- list(LRR_mean = list(CN_1 = -0.4156184, CN_3 = 0.1734862),
                  LRR_sd   = list(CN_1 = 0.2502591*sd_par, CN_3 = 0.2249798*sd_par))  # sd for one SNP

