#!/usr/bin/env Rscipt

args <- commandArgs( trailingOnly = TRUE )
wk_dir  <- args[1] ## path to IPQ.stats.txt generated in step (1)

suppressMessages({
  require(ggplot2)
  require(cowplot)
})

# PCA --------------------------------------------------------------------

dat <- readRDS(file = file.path(wk_dir, "IPQ.stats.txt"))

idx1 <- which( names(dat) == "Sample_ID" )
dat_pca <- dat[, -idx1]
mat <- as.matrix(dat_pca)
rownames <- dat$Sample_ID

PCA <- prcomp(mat, scale. = TRUE)
PC  <- predict(PCA)

PC  <- data.frame(Sample_ID = rownames(PC), 
                  PC, 
                  stringsAsFactors = FALSE)

write.table(PC, file = file.path(wk_dir, "IPQ_stats_PCA_res.txt"),
            quote = F, row.names = F, sep = "\t")

p12 <- ggplot() + 
  geom_point(data = PC, aes(PC1, PC2), shape = 1, size = 3) + 
  xlab("PC1") + 
  ylab("PC2") + 
  theme_bw(base_size = 9)+
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 20, hjust = 0.5)) + 
  ggtitle("PC2 ~ PC1")

p13 <- ggplot() + 
  geom_point(data = PC, aes(PC1, PC3), shape = 1, size = 3) + 
  xlab("PC1") + 
  ylab("PC3") + 
  theme_bw(base_size = 9)+
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 20, hjust = 0.5)) + 
  ggtitle("PC3 ~ PC1")

p23 <- ggplot() + 
  geom_point(data = PC, aes(PC2, PC3), shape = 1, size = 3) + 
  xlab("PC2") + 
  ylab("PC3") + 
  theme_bw(base_size = 9)+
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 20, hjust = 0.5)) + 
  ggtitle("PC3 ~ PC2")

png(filename = file.path(wk_dir, "IPQ_stats_PCA_plots.png"),
    width = 12, height = 12, units = "in", res = 512)
p <- plot_grid(p12, p13, p23, nrow = 2, labels = LETTERS[1:3])
print(p)
dev.off()






