rm(list=ls())
library(Seurat)
library(Signac)
library(tiltedCCA)

load("../../../../out/main/10x_greenleaf_timeseries.RData")

colnames(rna_predicted_mat) <- colnames(rna_common)
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
mentioned_genes <- mentioned_genes[mentioned_genes %in% colnames(rna_common)]

n <- nrow(rna_common)
for(k in 1:2){
  genes <- mentioned_genes[((k-1)*14+1):min((k*14), length(mentioned_genes))]
  
  png(paste0("../../../../out/figures/Writeup14o/Writeup14o_10x_greenleaf_timeseries_npregfast_mentioned_",
             k, ".png"),
      height = 2500, width = 3000, units = "px", res = 300)
  par(mfrow = c(3,5), mar = c(4,4,4,0.5))
  for(gene in genes){
    set.seed(10)
    x_vec <- c(1:n, 1:n)
    y_vec <- c(rna_common[,gene], atac_pred[,gene])
    col <- c(rep(rgb(0.5, 0.5, 0.5, 0.2), n), rep(rgb(0.8, 0, 0, 0.2), n))
    idx <- sample(1:length(x_vec))
    plot(x = x_vec[idx], y = y_vec[idx], 
         pch = 16, col = col[idx], main = gene)
    lines(rna_predicted_mat[,gene], col = "white", lwd = 4)
    lines(rna_predicted_mat[,gene], col = 1, lwd = 2)
    
    lines(atac_predicted_mat[,gene], col = "white", lwd = 4)
    lines(atac_predicted_mat[,gene], col = 2, lwd = 2)
  }
  graphics.off()
}

