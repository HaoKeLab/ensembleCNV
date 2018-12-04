#!/usr/bin/env Rscript

suppressMessages({
  require(data.table, quietly = TRUE)
  require(tibble, quietly = TRUE)
  require(optparse, quietly = TRUE)
})

option_list = list(
  make_option(c("-i", "--input"), action = "store", default = NA, type = "character",
              help = "path of perl code output"),
  make_option(c("-o", "--output"), action = "store", default = NA, type = "character",
              help = "path to save .rds file"),
  make_option(c("-s", "--startChr"), action = "store", default = 1, type = "integer",
              help = "start Chr name [default %default]"),
  make_option(c("-d", "--endChr"), action = "store", default = 22, type = "integer",
              help = "end Chr name [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))
pars <- c(opt$input, opt$output, opt$startChr, opt$endChr)

if ( any(is.na(pars)) ) {
  stop("All parameters must be supplied. (--help for detail)")
}

path_input <- opt$input
path_output <- opt$output
startChr <- opt$startChr
endChr <- opt$endChr

if ( !(is.integer(startChr) & is.integer(endChr)) ) {
  stop("parameters startChr and endChr must be integer.")
}

if ( startChr > endChr | startChr < 1 | endChr > 22 ) {
  stop("parameters startChr and endChr should satisfy 1 <= startChr < endChr <= 22.")
}

chrs <- seq(startChr, endChr)
# create LRR/BAF folder ---------------------------------------------------/
if ( !dir.exists(file.path(path_output, "LRR"))) dir.create(path = file.path(path_output, "LRR"), showWarnings = FALSE, recursive = TRUE)
if ( !dir.exists(file.path(path_output, "BAF"))) dir.create(path = file.path(path_output, "BAF"), showWarnings = FALSE, recursive = TRUE)

# read in annotate files --------------------------------------------------/
dat_snpName = fread( input = file.path(path_input, "snps_name.txt"), header = FALSE)
dat_snpName = as.data.frame(dat_snpName, stringsAsFactors = FALSE)
names( dat_snpName) <- c("Chr", "SNPs")

dat_snpNum = fread( input = file.path(path_input, "snps_number.txt"), header = FALSE)
dat_snpNum = as.data.frame( dat_snpNum, stringsAsFactors = FALSE)
names( dat_snpNum) <- c("Chr", "number")

dat_snpPos = fread( input = file.path(path_input, "SNP_pos.txt"), header = TRUE)
dat_snpPos = as.data.frame(dat_snpPos, stringsAsFactors = FALSE)
names( dat_snpPos) <- c("name", "chr","position")

dat_samples_order <- read.table(file = file.path(path_input, "samples_order.txt"), 
                                sep = "\t", header = F, stringsAsFactors = F)
names(dat_samples_order) <- c("sampleID" ,"order")
dat_samples_order <- dat_samples_order[order(dat_samples_order$order), ]

for (chr1 in chrs) {
  
  cat("chr:", chr1, "\n")
  
  snp1 <- unlist(strsplit( dat_snpName$SNPs[dat_snpName$Chr == chr1], 
                           split = "___", fixed = TRUE))
  
  n1 <- dat_snpNum$number[ dat_snpNum$Chr == chr1]
  stopifnot( length(snp1) == n1)
  
  snp_position_chr1 <- subset( dat_snpPos, name %in% snp1)
  snp_position_chr1 <- snp_position_chr1[ order(snp_position_chr1$position), ]
  
  snp1_order <- snp_position_chr1$name
  stopifnot( nrow(snp_position_chr1) == n1 )
  
  # read in LRR/BAF
  dat_chr1_LRR <- fread( input = file.path( path_input, "LRR", paste0(chr1, ".tab")), header = FALSE)
  dat_chr1_BAF <- fread( input = file.path( path_input, "BAF", paste0(chr1, ".tab")), header = FALSE)
  
  dat_chr1_LRR <- as.data.frame(dat_chr1_LRR, stringsAsFactors = FALSE)
  dat_chr1_BAF <- as.data.frame(dat_chr1_BAF, stringsAsFactors = FALSE)
  
  rownames( dat_chr1_LRR ) <- NULL
  rownames( dat_chr1_BAF ) <- NULL
  
  dat_chr1_LRR <- column_to_rownames( dat_chr1_LRR, var = "V1")
  dat_chr1_BAF <- column_to_rownames( dat_chr1_BAF, var = "V1")
  
  stopifnot( ncol(dat_chr1_LRR) == n1 )
  stopifnot( ncol(dat_chr1_BAF) == n1 )
  
  names(dat_chr1_LRR) <- snp1
  names(dat_chr1_BAF) <- snp1
  
  dat_chr1_LRR <- dat_chr1_LRR[, snp1_order, drop = FALSE]
  dat_chr1_BAF <- dat_chr1_BAF[, snp1_order, drop = FALSE]
  
  ## check samples_order
  stopifnot( all(rownames(dat_chr1_LRR) == dat_samples_order$sampleID) )
  stopifnot( all(rownames(dat_chr1_BAF) == dat_samples_order$sampleID) )
  
  saveRDS( dat_chr1_LRR, file = file.path(path_output, "LRR", paste0("matrix_chr_", chr1, "_LRR.rds")))
  saveRDS( dat_chr1_BAF, file = file.path(path_output, "BAF", paste0("matrix_chr_", chr1, "_BAF.rds")))
  
}

cat("Analysis completed! The output files are at:", path_output, "\n")

