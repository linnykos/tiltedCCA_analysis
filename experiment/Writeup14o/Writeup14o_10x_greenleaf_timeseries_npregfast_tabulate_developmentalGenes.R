rm(list=ls())
library(Seurat)
library(Signac)
library(tiltedCCA)

load("../../../../out/main/10x_greenleaf_timeseries.RData")
load("../../../../out/main/10x_greenleaf_developmentalGenes.RData")

colnames(rna_predicted_mat) <- colnames(rna_common)
colnames(atac_predicted_mat) <- colnames(atac_pred)

n <- nrow(rna_common)
for(k in 1:2){
  genes <- selection_res$selected_variables[((k-1)*25+1):min((k*25), length(selection_res$selected_variables))]
  
  png(paste0("../../../../out/figures/Writeup14o/Writeup14o_10x_greenleaf_timeseries_npregfast_developmental_",
             k, ".png"),
      height = 3500, width = 3000, units = "px", res = 300)
  par(mfrow = c(5,5), mar = c(4,4,4,0.5))
  for(gene in genes){
    set.seed(10)
    x_vec <- c(1:n, 1:n)
    y_vec <- c(rna_common[,gene], atac_pred[,gene])
    col <- c(rep(rgb(0.5, 0.5, 0.5, 0.2), n), rep(rgb(0.8, 0, 0, 0.2), n))
    idx <- sample(1:length(vec))
    plot(x = x_vec[idx], y = y_vec[idx], 
         pch = 16, col = col[idx], main = gene)
    lines(rna_predicted_mat[,gene], col = "white", lwd = 4)
    lines(rna_predicted_mat[,gene], col = 1, lwd = 2)
    
    lines(atac_predicted_mat[,gene], col = "white", lwd = 4)
    lines(atac_predicted_mat[,gene], col = 2, lwd = 2)
  }
  graphics.off()
}

