

path_output <- ""
file_ipattern <- ""
file_penncnv <- ""
file_quantisnp <- ""

col_sel <- c("chr", "posStart", "posEnd", "CN", "Sample_ID", "conf", 
             "numSNP", "avgConf", "length", "CNV_type", "method") ## add method

# generate data from iPattern, PennCNV, QuantiSNP  ------------------------

## ipattern ---------------------------------------------------------------------

read_icnv <- function(file_icnv, col_sel) {
  
  library(data.table)
  cat("Read in CNV from iPattern.\n")
  dat <- read.table(file = file_icnv, sep = "\t", skip = 17, check.names = FALSE, 
                    header = FALSE, comment.char = "",
                    stringsAsFactors = FALSE) # check if snpnames have "#"
  
  names(dat) <- c("CNV_type", "chr", "posStart", "posEnd",
                  "numSNP", "on_probe.num", "clusterIdx",
                  "gain_loss_score", "cluster_score", "gain_loss_sample.num",
                  "conf", "Sample_ID", "CNV_event_ID", "CNVR_ID")
  dat$length <- dat$posEnd - dat$posStart + 1
  dat$avgConf <- dat$conf/dat$numSNP
  
  dat$Sample_ID <- gsub("^X", "", dat$Sample_ID, perl = TRUE)
  dat$Sample_ID <- gsub(".", "-", dat$Sample_ID, fixed = TRUE)
  
  # filter chr, CNV_type
  dat <- subset(dat, chr %in% c(1:22) & CNV_type %in% c("Gain", "Loss"))
  # dat$Sample_ID <- gsub(pattern = "X", replacement = "", dat$Sample_ID, fixed = TRUE)  ## check
  dat$CN <- ifelse(dat$CNV_type == "Gain", 3, 1) ##
  dat$chr <- as.integer(dat$chr)
  dat$method <- "iPattern"
  
  dat[, col_sel]  ## select columns
}

# merge all groups results
dat_ipattern <- read_icnv( file_icnv = file_ipattern, col_sel = col_sel )

saveRDS( dat_ipattern, file = file.path(path_output, "cnvs.ipattern.rds"))


# penncnv -----------------------------------------------------------------

read_pcnv <- function(file_pcnv, col_sel) {
  
  cat("Read in CNV from PennCNV.\n")
  dat <- read.table(file = file_pcnv, sep = "\t", check.names = FALSE,
                    header = FALSE, stringsAsFactors = FALSE, 
                    comment.char = "")
  names(dat) <- c("chr", "posStart", "posEnd", "CN", "Sample_ID", "snpStart", "snpEnd", "conf", "numSNP")
  dat$Sample_ID <- gsub(".txt", "", dat$Sample_ID)
  dat$length <- dat$posEnd - dat$posStart + 1
  dat$avgConf <- dat$conf/dat$numSNP
  dat$CNV_type <- ifelse(dat$CN > 2, "Gain", "Loss")
  dat$method <- "PennCNV"
  dat$CN[which(dat$CN >= 3)] <- 3 # set CN >= 3 all = 3
  
  dat[, col_sel]
}

dat_penncnv <- read_pcnv(file_pcnv = file_penncnv, col_sel = col_sel)

saveRDS( dat_penncnv, file = file.path(path_output, "cnvs.penncnv.rds"))


# quantisnp ---------------------------------------------------------------

# generate quantisnp.cnv use perl code

read_qcnv <- function(file_qcnv, col_sel) {
  
  cat("Read in CNV from QuantiSNP.\n")
  dat <- read.table(file = file_qcnv, 
                    sep = "\t", 
                    header = TRUE, 
                    check.names = FALSE, 
                    stringsAsFactors = FALSE,
                    comment.char = "")
  ## Max.log BF
  names(dat) <- c("Sample_ID", "chr", "posStart", "posEnd", "snpStart", "snpEnd", "length", "numSNP", 
                  "CN", "conf", "Log_BF.State.0", "Log_BF.State.1", "Log_BF.State.2", "Log_BF.State.3",
                  "Log_BF.State.4", "Log_BF.State.5", "Log_BF.State.6")
  dat <- subset(dat, chr %in% c(1:22))
  dat$CNV_type <- ifelse(dat$CN > 2, "Gain", "Loss")
  dat$avgConf <- dat$conf/dat$numSNP
  dat$method <- "QuantiSNP"
  dat$CN[which(dat$CN >= 3)] <- 3  # set CN >= 3 all = 3
  
  dat[, col_sel]
}

dat_quantisnp <- read_qcnv(file_qcnv = file_quantisnp, col_sel = col_sel)

saveRDS( dat_quantisnp, file = file.path(path_output, "cnvs.quantisnp.rds"))

