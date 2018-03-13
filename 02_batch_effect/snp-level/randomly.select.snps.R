

## test on starnet
require(data.table)

path_snps <- "/sc/orga/projects/haok01a/chengh04/Food_Allergy/raw_data_batch/batch1"
dat_snps <- fread(input = file.path(path_snps, "FA_batch1_recluster_SNP_Table.txt"))

path_output <- "/sc/orga/projects/haok01a/chengh04/Food_Allergy/code_batch/SNP_recluster/dat"

dat_snps <- as.data.frame(dat_snps, stringsAsFactors = FALSE)
dat_snps <- subset(dat_snps, Chr %in% 1:22)
set.seed(123456789)
snps <- sample( dat_snps$Name )
snps <- snps[1:100000]

write.table(snps, file = file.path(path_output, "snps.randomly.select.txt"),
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
