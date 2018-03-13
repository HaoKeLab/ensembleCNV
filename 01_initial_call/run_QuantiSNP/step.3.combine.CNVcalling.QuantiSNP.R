#!/hpc/packages/minerva-common/R/3.3.1/lib64/R/bin/Rscript --vanilla

suppressPackageStartupMessages(require(optparse))
require(data.table)

option_list <- list(
  make_option(c("-r", "--res"), action = "store", default = NA, type = "character",
              help = "folder of QuantiSNP results."),
  make_option(c("-o", "--output"), action = "store", default = NA, type = "character",
              help = "output the combined result."),
  make_option(c("-n", "--name"), action = "store", default = NA, type = "character",
              help = "project name")
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.na(opt$res) | is.na(opt$output) | is.na(opt$name)) {
  stop("all three args must be supplied.")
}

path_res <- opt$res
path_output <- opt$output
name_projet <- opt$name

## combine all samples cnv file for all columns
# for 2017
# path_res <- "/sc/orga/projects/haok01a/chengh04/FMD.GWAS/callCNV/Kovacic_128samples_042717/QuantiSNP/res"
# path_output <- "/sc/orga/projects/haok01a/chengh04/FMD.GWAS/callCNV/Kovacic_128samples_042717/QuantiSNP"
# name_projet <- "Kovacic_128samples_042717"

# for 2015
# path_res <- "/sc/orga/projects/haok01a/chengh04/FMD.GWAS/callCNV/Valentina_112samples_120214/QuantiSNP/res"
# path_output <- "/sc/orga/projects/haok01a/chengh04/FMD.GWAS/callCNV/Valentina_112samples_120214/QuantiSNP"
# name_projet <- "Valentina_112samples_120214"

folders  <- list.files(path = path_res)

res <- data.frame()
for (i in 1:length(folders)) {
  
  folder1 <- folders[i]
  path1   <- file.path(path_res, folder1)
  files <- list.files(path = path1)
  
  cat("Read in:", path1, "\n")
  dat1 <- fread(input = file.path(path1, paste0(folder1, ".cnv")), check.names = FALSE)
  dat1 <- as.data.frame(dat1, stringsAsFactors = FALSE)
  
  res <- rbind(res, dat1)
}

write.table(res, file = file.path(path_output, paste0(name_projet, ".cnv")),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
