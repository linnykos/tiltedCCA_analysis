rm(list=ls())
library(Seurat)
library(Signac)
library(tiltedCCA)

load("../../../../out/main/10x_greenleaf_timeseries.RData")
load("../../../../out/main/10x_greenleaf_developmentalGenes.RData")

colnames(rna_predicted_mat) <- colnames(rna_common)
colnames(atac_predicted_mat) <- colnames(atac_pred)

n <- nrow(rna_common)
base_palette <- RColorBrewer::brewer.pal(11, name = "RdYlBu")
color_vec <- rev(grDevices::colorRampPalette(base_palette)(n))

rna_common <- t(apply(rna_common, 1, function(x){scale(x)}))
colnames(rna_common) <- colnames(rna_predicted_mat)
atac_pred <- t(apply(atac_pred, 1, function(x){scale(x)}))
colnames(atac_pred) <- colnames(atac_predicted_mat)
rna_predicted_mat <- t(apply(rna_predicted_mat, 1, function(x){scale(x)}))
colnames(rna_predicted_mat) <- colnames(rna_common)
atac_predicted_mat <- t(apply(atac_predicted_mat, 1, function(x){scale(x)}))
colnames(atac_predicted_mat) <- colnames(atac_pred)

mentioned_genes <- c("ALDH2", "APOE", "AQP4", "ASCL1", "BHLHE22",
                     "BHLHE40", "C16orf89", "CAV2", "DLX2", "DOK5",
                     "DUSP1", "EOMES", "ETV4", "FOS", "FOXJ1", "GLI3",
                     "HAS2", "HES1", "HES4", "HSPA1A", "HSPA1B",
                     "ID3", "IGFBP7", "JUN", "KIF1A", "LIMCH1",
                     "MBP", "MEF2C", "NEUROD1", "NEUROD2", "NEUROD4",
                     "NEUROD6", "NEUROG1", "NEUROG2", "NFIA", "NFIB", "NFIC",
                     "NHLH1", "NR2F1", "PAX6", "RFX4", "RUNX1",
                     "OLIG1", "OLIG2", "SOX2", "SOX3", "SOX6",
                     "SOX9", "SOX10", "SOX21", "SPARCL1", "SNCB", "TBX",
                     "TNC", "TOP2A", "TRB1", "WNT11")
mentioned_genes <- intersect(mentioned_genes, colnames(rna_common))

genes <- mentioned_genes 
png(paste0("../../../../out/figures/Writeup14p/Writeup14p_10x_greenleaf_leafplots-developmental_1.png"),
    height = 3500, width = 3000, units = "px", res = 300)
par(mfrow = c(6,5), mar = c(4,4,4,0.5))
for(gene in genes){
  set.seed(10)
  idx <- sample(1:n)
  
  x_vec1 <- atac_pred[idx,gene]
  x_vec2 <- atac_predicted_mat[idx,gene]
  y_vec1 <- rna_common[idx,gene]
  y_vec2 <- rna_predicted_mat[idx,gene]
  
  x_mean <- mean(x_vec2); x_vec1 <- x_vec1- x_mean; x_vec2 <- x_vec2 - x_mean
  x_sd <- sd(x_vec2); x_vec1 <- x_vec1/x_sd; x_vec2 <- x_vec2/x_sd
  y_mean <- mean(y_vec2); y_vec1 <- y_vec1- y_mean; y_vec2 <- y_vec2 - y_mean
  y_sd <- sd(y_vec2); y_vec1 <- y_vec1/y_sd; y_vec2 <- y_vec2/y_sd
  
  plot(x = x_vec1, y = y_vec1,
       pch = 16, col = color_vec[idx], main = gene, cex = 0.5,
       xlab = "ATAC", ylab = "RNA", asp = T)
  
  points(x = x_vec2, 
         y = y_vec2, 
         pch = 16, col = "white", cex = 2)
  points(x = x_vec2, 
         y = y_vec2, 
         pch = 16, col = color_vec[idx], cex = 1.5)
  
  points(x = x_vec2[1], 
         y = y_vec2[1], 
         pch = 16, col = "white", cex = 3)
  points(x = x_vec2[1], 
         y = y_vec2[1], 
         pch = 16, col = color_vec[1], cex = 2.5)
  
  x_right <- mean(x_vec1[x_vec1 >= quantile(x_vec1, probs = 0.9)])
  x_left <- mean(x_vec1[x_vec1 <= quantile(x_vec1, probs = 0.1)])
  y_top <- mean(y_vec1[y_vec1 >= quantile(y_vec1, probs = 0.9)])
  y_bot <- mean(y_vec1[y_vec1 <= quantile(y_vec1, probs = 0.1)])
  lines(c(2*x_right,2*x_left), c(2*y_top, 2*y_bot), lwd = 2, lty = 2)
}
graphics.off()

############################################################

