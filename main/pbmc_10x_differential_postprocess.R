rm(list=ls())
load("../../../out/main/10x_pbmc_differential.RData")
load("../../../out/main/10x_pbmc_tiltedcca.RData")

library(Seurat)
library(Signac)
library(tiltedCCA)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

Seurat::DefaultAssay(pbmc) <- "SCT"
logpval_vec <- sapply(Seurat::VariableFeatures(pbmc), function(gene){
  val <- stats::median(
    sapply(1:length(gene_de_list), function(j){
      idx <- which(rownames(gene_de_list[[j]]) == gene)
      if(length(idx) == 0) return(1)
      gene_de_list[[j]][idx, "p_val"]
    })
  )
  
  -log10(val)
})
names(logpval_vec) <- Seurat::VariableFeatures(pbmc)

rsquare_vec <- tiltedCCA:::postprocess_alignment(multiSVD_obj = multiSVD_obj,
                                                 bool_use_denoised = T,
                                                 seurat_obj = pbmc,
                                                 multiSVD_assay = 1,
                                                 seurat_assay = "SCT",
                                                 seurat_slot = "data")
all(names(logpval_vec) == names(rsquare_vec))
stats::median(rsquare_vec[which(logpval_vec >=1)])

png("../../../out/figures/main/10x_pbmc_differential_gene.png",
    height = 3500, width = 2500, res = 500, units = "px")
par(mar = c(5,5,4,1))
tiltedCCA:::plot_alignment(rsquare_vec = rsquare_vec,
                           logpval_vec = logpval_vec,
                           main = "PBMC (10x, RNA+ATAC)\nGene differentiability vs. alignment",
                           bool_mark_ymedian = T,
                           col_gene_highlight_border = rgb(255, 205, 87, 255*0.5, maxColorValue = 255),
                           col_points = rgb(0.5, 0.5, 0.5, 0.1),
                           cex_axis = 1.5, 
                           cex_lab = 1.5,
                           cex_points = 2.5,
                           lty_polygon = 2,
                           lwd_grid_major = 2,
                           lwd_grid_minor = 1,
                           lwd_axis = 1.5,
                           lwd_axis_ticks = 1.5,
                           lwd_polygon_bold = 5)
graphics.off()
