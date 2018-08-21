#!/usr/bin/env Rscript

# generate the list of 100000 randomly selected SNPs
# Input is SNP_Map.txt generated from Genome Studio
# including at least 3 columns: Name, Chromosome, Position

suppressMessages({
  require( data.table, quietly = TRUE)
})

args <- commandArgs( trailingOnly = TRUE )
file_snps <- args[1]   ## SNP_Table.txt file from Genome Studio
path_output <- args[2] ## path to save randomly selected SNPs

## sampleing from chr: 1-22
dat_snps <- fread( input = file_snps )
dat_snps <- as.data.frame(dat_snps, stringsAsFactors = FALSE)
dat_snps <- subset(dat_snps, Chromosome %in% 1:22)

snps <- sample( dat_snps$Name )
snps.selected <- snps[ 1:100000 ]

write.table(snps.selected, file = file.path(path_output, "snps.down.sample.txt"),
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
