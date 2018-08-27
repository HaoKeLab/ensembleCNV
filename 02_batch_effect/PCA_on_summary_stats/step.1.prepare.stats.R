#!/urs/bin/env Rscript

args <- commandArgs( trailingOnly = TRUE )

path_ipattern  <- args[1]
path_penncnv   <- args[2]
path_quantisnp <- args[3]
path_output    <- args[4]

suppressMessages({
  require(data.table)
})

# ipattern ----------------------------------------------------------------
## for number of samples larger than 500, samples may need to be splited into batches to run ipattern
read_ipattern_batch <- function(path_ipattern) {
  
  ## NumCNV
  cnv_file <- list.files(path = path_ipattern, pattern = "_all_calls.txt$")
  
  dat <- read.table(file = file.path(path_ipattern, cnv_file),
                    header = FALSE, sep = "\t", comment.char = "#", as.is = TRUE)
  names(dat) <- c("CNV_type", "chr", "posStart", "posEnd",
                  "numSNP", "on_probe.num", "clusterIdx",
                  "gain_loss_score", "cluster_score", "gain_loss_sample.num",
                  "conf", "Sample_ID", "CNV_event_ID", "CNVR_ID")
  dat <- subset(dat, chr %in% c(1:22))
  tbl <- table(dat$Sample_ID)
  dat_tbl <- as.data.frame(tbl)
  names(dat_tbl) <- c("Sample_ID", "iPattern.NumCNV")
  NumSample <- nrow(dat_tbl)
  
  ## sample.stats.txt
  stat_file <- list.files(path = path_ipattern, pattern = "_sample.stats.txt$")
  
  dat_stat <- read.table(file = file.path(path_ipattern, stat_file),
                         header = FALSE, sep = "\t", nrows = NumSample, as.is = TRUE)
  names(dat_stat) <- c("Sample_ID", "iPattern.LRR_SD", "iPattern.base_CN")
  dat_stat <- dat_stat[, c("Sample_ID", "iPattern.LRR_SD")]
  
  ## clean sample ID: remove path information, remove subfix ".rescale"
  samples <- dat_stat$Sample_ID
  samples <- unlist( lapply(1:length(samples), FUN = function(k) {
    sample1 <- samples[k]
    strs <- unlist(strsplit(sample1, split = "/", fixed = TRUE))
    str1 <- strs[length(strs)]
  }) )
  samples <- gsub("\\.rescale$", "", samples)
  dat_stat$Sample_ID <- samples
  
  res <- merge(dat_stat, dat_tbl)
  ## if Sample_ID starts with number
  res$Sample_ID <- gsub(pattern = "^X", replacement = "", res$Sample_ID, perl = TRUE)  ## check
  
  res
}

cat("Processing iPattern results ...\n")
dat_stats_ipattern <- read_ipattern_batch(path_ipattern = path_ipattern)

write.table(dat_stats_ipattern, 
            file = file.path(path_output, "ipattern.stats.txt"),
            quote = F, row.names = F, sep = "\t")
cat("Done.\n")

# penncnv sample-level -----------------------------------------------------
cat("Processing PennCNV results ...\n")
dat_penncnv <- read.table(file = file.path(path_penncnv, "CNV.PennCNV_qc_new.txt"),
                          sep = "\t",
                          header = TRUE,
                          check.names = FALSE,
                          stringsAsFactors = FALSE)
dat_penncnv$File <- gsub("\\.txt$", "", dat_penncnv$File, perl = TRUE) 
dat_penncnv$WF <- abs(dat_penncnv$WF)

fp <- c( "LRR_SD", "BAF_SD", "BAF_drift", "WF", "NumCNV" )
dat_penncnv <- dat_penncnv[, c("File", fp)]
names(dat_penncnv) <- c("Sample_ID", paste("PennCNV", fp, sep = "."))

dat_stats_penncnv <- dat_penncnv

write.table(dat_stats_penncnv, 
            file = file.path(path_output, "penncnv.stats.txt"),
            quote = F, row.names = F, sep = "\t")
cat("Done.\n")

# quantisnp ---------------------------------------------------------------
read_quantisnp_per_sample <- function(path_res, sample_id) {
  
  ## get numCNV
  file_cnv <- file.path(file.path(path_res, sample_id), paste0(sample_id, ".cnv"))
  dat_cnv <- fread(input = file_cnv)
  numCNV <- sum(dat_cnv$Chromosome %in% c(1:22))
  
  ## get LRR.SD and BAF.SD
  ## Note: in the .qc file, QuantiSNP has formatting issue
  ##       the column name "Gender" is written at the start of the second line
  file_qc <- file.path(file.path(path_res, sample_id), paste0(sample_id, ".qc"))
  dat_line2 <- read.table(file = file_qc, skip = 1, nrows = 1, header = FALSE, stringsAsFactors = FALSE)
  dat_line2 <- dat_line2[, -1]
  names(dat_line2) <- c("Sample_ID", "Chr", "OutlierRate", "LRR_SD", "BAF_SD", "Gender")
  
  dat_other <- read.table(file = file_qc, skip = 2, header = FALSE, stringsAsFactors = FALSE)
  names(dat_other) <- c("Sample_ID", "Chr", "OutlierRate", "LRR_SD", "BAF_SD", "Gender")
  dat <- rbind(dat_line2, dat_other)
  names(dat) <- c("Sample_ID", "Chr", "OutlierRate", "LRR_SD", "BAF_SD", "Gender")
  dat <- subset(dat, Chr %in% c(1:22))
  
  Sample_ID <- unique(dat$Sample_ID)
  LRR_SD <- mean(dat$LRR_SD, na.rm = TRUE)
  BAF_SD <- mean(dat$BAF_SD, na.rm = TRUE)
  
  res1 <- data.frame(Sample_ID = Sample_ID, 
                     QuantiSNP.NumCNV = numCNV,
                     QuantiSNP.LRR_SD = LRR_SD, 
                     QuantiSNP.BAF_SD = BAF_SD,
                     stringsAsFactors = FALSE)
  return(res1) ## for one sample
}

read_quantisnp <- function(path_res) {
  
  samples <- list.files(path = path_res)
  res <- data.frame() ## all QuantiSNP statistics
  for (i in 1:length(samples)) {
    
    sample1 <- samples[i]
    #cat("i:", i, length(samples), "SampleID:", sample1, "\n")
    
    res1 <- read_quantisnp_per_sample(path_res = path_res, sample_id = sample1)
    res <- rbind(res, res1)
  }
  res
}

cat("Processing QuantiSNP results ...\n")
dat_stats_quantisnp <- read_quantisnp(path_res = path_quantisnp)

write.table(dat_stats_quantisnp, 
            file = file.path(path_output, "quantisnp.stats.txt"),
            quote = F, row.names = F, sep = "\t")
cat("Done.")


# IPQ ---------------------------------------------------------------------

cat("Combine summary statistics from different methods ...\n")
## iPattern converts "-" in Sample_ID to "."
## recover the original Sample_ID
idx <- grep("-", dat_stats_penncnv$Sample_ID) 
samples.raw <- dat_stats_penncnv$Sample_ID[ idx  ]
samples.alt <- sub("-", ".", samples.raw)

for (i in 1:length(samples.alt)) {
	idx1 <- which(dat_stats_ipattern$Sample_ID == samples.alt[i])
	dat_stats_ipattern$Sample_ID[ idx1 ] <- samples.raw[i]
}

res_IP <- merge(dat_stats_ipattern, dat_stats_penncnv)
stopifnot( nrow(res_IP) == nrow(dat_stats_ipattern))

res_IPQ <- merge(res_IP, dat_stats_quantisnp)
stopifnot( nrow(res_IPQ) == nrow(res_IP) )

write.table(res_IPQ, 
            file = file.path(path_output, "IPQ.stats.txt"),
            quote = F, row.names = F, sep = "\t")
cat("Done.\n")




