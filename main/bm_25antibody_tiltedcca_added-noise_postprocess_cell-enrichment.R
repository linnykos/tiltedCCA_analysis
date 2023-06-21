rm(list=ls())
library(Seurat)
library(Signac)
library(tiltedCCA)
source("bm_25antibody_colorPalette.R")

percentage_vec <- c(0, 0.1, 0.25, 0.5, 1) # c(0.1, 0.25, 0.5, 1)
for(percentage in percentage_vec){
  print(paste0("Trying percentage ", percentage))
  load(paste0("../../../out/main/citeseq_bm25_tcca_added-noise_", percentage, ".RData"))
  
  print("Computing")
  membership_vec <- factor(bm$celltype.l2)
  cell_enrichment_res <- tiltedCCA:::postprocess_cell_enrichment(
    input_obj = multiSVD_obj, 
    membership_vec = membership_vec, 
    max_subsample = 1000,
    verbose = 1
  )
  
  y_vec <- cell_enrichment_res$enrichment_common$df[,"value"]
  names(y_vec) <- cell_enrichment_res$enrichment_common$df[,"celltype"]
  x_vec <- cell_enrichment_res$enrichment_distinct_2$df[,"value"]
  names(x_vec) <- cell_enrichment_res$enrichment_distinct_2$df[,"celltype"]
  
  png(paste0("../../../out/figures/main/bm_25antibody_added-noise_", percentage, "_cell_enrichment_common-vs-distinct2.png"),
      height = 2500, width = 2500, res = 500, units = "px")
  par(mar = c(4, 4, 0.5, 0.5))
  plot(NA, xlim = c(0,1), ylim = c(0,1),
       main = "", xlab = "", ylab = "",
       xaxt = "n", yaxt = "n", bty = "n", asp = T)
  for(x in seq(0,1,by=0.1)){
    lines(rep(x,2), c(-10,10), col = "gray", lty = 3, lwd = 2)
  }
  for(y in seq(0,1,by=0.1)){
    lines(c(-10,10), rep(y,2), col = "gray", lty = 3, lwd = 2)
  }
  lines(c(-10,10), c(-10,10), col = 2, lty = 2, lwd = 3)
  
  col_vec <- col_palette[names(x_vec)]
  points(x_vec, y_vec, col = 1, cex = 4, pch = 16)
  points(x_vec, y_vec, col = "white", cex = 3, pch = 16)
  points(x_vec, y_vec, col = col_vec, cex = 2.5, pch = 16)
  
  axis(1, cex.axis = 1.25, cex.lab = 1.25)
  axis(2, cex.axis = 1.25, cex.lab = 1.25)
  
  graphics.off()
  
  
  y_vec <- cell_enrichment_res$enrichment_common$df[,"value"]
  names(y_vec) <- cell_enrichment_res$enrichment_common$df[,"celltype"]
  x_vec <- -cell_enrichment_res$enrichment_distinct_1$df[,"value"]
  names(x_vec) <- cell_enrichment_res$enrichment_distinct_1$df[,"celltype"]
  
  png(paste0("../../../out/figures/main/bm_25antibody_added-noise_", percentage, "_cell_enrichment_common-vs-distinct1.png"),
      height = 2500, width = 2500, res = 500, units = "px")
  par(mar = c(4, 4, 0.5, 0.5))
  plot(NA, xlim = c(-1,0), ylim = c(0,1),
       main = "", xlab = "", ylab = "",
       xaxt = "n", yaxt = "n", bty = "n", asp = T)
  for(x in seq(-1,0,by=0.1)){
    lines(rep(x,2), c(-10,10), col = "gray", lty = 3, lwd = 2)
  }
  for(y in seq(0,1,by=0.1)){
    lines(c(-10,10), rep(y,2), col = "gray", lty = 3, lwd = 2)
  }
  lines(c(-10,10), c(10,-10), col = 2, lty = 2, lwd = 3)
  
  col_vec <- col_palette[names(x_vec)]
  points(x_vec, y_vec, col = 1, cex = 4, pch = 16)
  points(x_vec, y_vec, col = "white", cex = 3, pch = 16)
  points(x_vec, y_vec, col = col_vec, cex = 2.5, pch = 16)
  
  axis(1, cex.axis = 1.25, cex.lab = 1.25)
  # axis(2, cex.axis = 1.25, cex.lab = 1.25)
  graphics.off()
}
