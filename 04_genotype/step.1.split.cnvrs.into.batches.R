#!/hpc/packages/minerva-common/R/3.3.1/lib64/R/bin/Rscript --vanilla

## split CNVR in each chr into batches for submiting jobs 

suppressMessages(require(optparse))

option_list = list(
  make_option(c("-i", "--input"), action = "store", type = "character", default = NA,
              help = "CNVR dataset input"),
  make_option(c("-o", "--output"), action = "store", type = "character", default = NA,
              help = "CNVR dataset output"),
  make_option(c("-c", "--cnv"), action = "store", type = "character", default = NA,
              help = "CNV after cleaning dataset"),
  make_option(c("-n", "--num"), action = "store", type = "integer", default = NA,
              help = "number for each batch in each chromosome")
)

opt = parse_args( OptionParser(option_list = option_list) )

pars = c(opt$input, opt$output, opt$num, opt$cnv)
if ( any(is.na(pars)) ) {
  stop("All parameter must be supplied.(--help for detail)")
}


# main  -------------------------------------------------------------------

dt_cnvr = readRDS( file = opt$input )
n_cnvr  = nrow(dt_cnvr)

dt_cnvr = dt_cnvr[order(dt_cnvr$chr, dt_cnvr$posStart, dt_cnvr$posEnd), ]

cat('total cnvr number:', n_cnvr, "\n")

number_each_batch = as.integer( opt$num )  ## 200 default

# add raw Freq information 
dt_cnv = readRDS(file = opt$cnv)
nrow(dt_cnv)
tbl <- table(dt_cnv$CNVR_ID)
freqs <- as.vector(tbl)
dt_freq <- data.frame(CNVR_ID = names(tbl), Freq = freqs, stringsAsFactors = FALSE)

dt_cnvr <- merge(dt_cnvr, dt_freq, by = "CNVR_ID")
stopifnot( nrow(dt_cnvr) == n_cnvr)

# split batches in each chr
chrs <- sort(unique(dt_cnvr$chr))

dt_cnvr_new <- data.frame()
for (chr1 in chrs) {
  
  dt_cnvr1 <- subset(dt_cnvr, chr == chr1)
  idxs_batch <- 1:nrow(dt_cnvr1)
  
  n1 <- nrow(dt_cnvr1)
  n2 <- ceiling(n1/number_each_batch) 
  
  cat("chr:", chr1, "number of cnvrs:", n1, "\n")
  if (n2 == 1) {
    
    dt_cnvr1$batch <- 1
    dt_cnvr_new <- rbind(dt_cnvr_new, dt_cnvr1)
    
  } else {
    
    cuts <- cut(idxs_batch, breaks = n2, include.lowest = TRUE)
    cuts_index <- as.integer(cuts)
    dt_cnvr1$batch <- cuts_index
    
    dt_cnvr_new <- rbind(dt_cnvr_new, dt_cnvr1)
  }
  
}

saveRDS(dt_cnvr_new, file = opt$output)


