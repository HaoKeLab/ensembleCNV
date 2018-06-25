#!/usr/bin/env Rscipt

suppressMessages({
  require(ggplot2)
  require(cowplot)
})

path_input  <- ""
path_output <- ""
# PCA -------------------------------------------------------------------/

dat <- readRDS(file = file.path(path_input, "IPQ.sample.level.statics.rds"))

idx1 <- which( names(dat) == "Sample_ID" )
dat_pca <- dat[, -idx1]
mat = as.matrix(dat_pca)


PCA = prcomp(mat, scale. = TRUE)
PC  = predict(PCA)
PC  = as.data.frame(PC, stringsAsFactors = FALSE)

rownames(PC) <- dat$Sample_ID

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

png(filename = file.path(path_output, "PCA.IPQ.sample.level.statics.png"),
    width = 12, height = 12, units = "in", res = 512)
p <- plot_grid(p12, p13, p23, nrow = 2, labels = LETTERS[1:3])
print(p)
dev.off()






