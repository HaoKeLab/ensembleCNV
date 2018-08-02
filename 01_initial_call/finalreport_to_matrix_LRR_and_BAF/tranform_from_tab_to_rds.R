
args <- commandArgs(trailingOnly = TRUE)

## after using perl script to split finalreport into LRR and BAF 'tab' format data
# transform from .tab format to .rds


suppressMessages({
  require(data.table, quietly = TRUE)
  require(tibble, quietly = TRUE)
})

path_input <- args[1] ## path_output from perl code
path_output <- args[2] ## path to save .rds data
chr_start <- as.integer(args[3]) ## start chr
chr_end   <- as.integer(args[4]) ## end chr


# create LRR/BAF folder ---------------------------------------------------/
if ( !dir.exists(file.path(path_output, "LRR"))) dir.create(path = file.path(path_output, "LRR"))
if ( !dir.exists(file.path(path_output, "BAF"))) dir.create(path = file.path(path_output, "BAF"))

# read in annotate files --------------------------------------------------/
dat_snpName = fread( input = file.path(path_input, "snps_name.txt"), header = FALSE)
dat_snpName = as.data.frame(dat_snpName, stringsAsFactors = FALSE)
names( dat_snpName) <- c("Chr", "SNPs")

dat_snpNum = fread( input = file.path(path_input, "snps_number.txt"), header = FALSE)
dat_snpNum = as.data.frame( dat_snpNum, stringsAsFactors = FALSE)
names( dat_snpNum) <- c("Chr", "number")

dat_snpPos = fread( input = file.path(path_input, "snps_position.txt"), header = FALSE)
dat_snpPos = as.data.frame(dat_snpPos, stringsAsFactors = FALSE)
names( dat_snpPos) <- c("name", "position")

for ( chr1 in chr_start:chr_end ) {
  
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
  
  names(dat_chr1_LRR) <- snp1
  names(dat_chr1_BAF) <- snp1
  
  dat_chr1_LRR <- dat_chr1_LRR[, snp1_order]
  dat_chr1_BAF <- dat_chr1_BAF[, snp1_order]
 
  saveRDS( dat_chr1_LRR, file = file.path(path_output, "LRR", paste0("matrix_chr_", chr1, "_LRR.rds")))
  saveRDS( dat_chr1_BAF, file = file.path(path_output, "BAF", paste0("matrix_chr_", chr1, "_BAF.rds")))
  
}


















