rm(list=ls())
load("../../../out/main/10x_mouseembryo_differential.RData")
load("../../../out/main/10x_mouseembryo_tiltedcca.RData")

library(Seurat)
library(Signac)
library(tiltedCCA)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

# temporary fix
num_celltypes <- ceiling(sqrt(2*length(gene_de_list)))
combn_mat <- combn(num_celltypes, 2)
stopifnot(ncol(combn_mat) == length(gene_de_list))
gene_de_list <- list(de_list = gene_de_list, combn_mat = combn_mat, level_vec = 1:num_celltypes)
############

Seurat::DefaultAssay(mbrain) <- "SCT"
gene_names <- Seurat::VariableFeatures(mbrain)
logpval_vec <- sapply(1:length(gene_names), function(k){
  if(k %% floor(length(gene_names)/10) == 0) cat('*')
  gene <- gene_names[k]
  
  # cycle through all the celltypes
  celltype_vec <- sapply(1:length(gene_de_list$level_vec), function(i){
    idx <- which(gene_de_list$combn_mat == i, arr.ind = T)[,2]
    vec <-  sapply(idx, function(j){
      idx <- which(rownames(gene_de_list$de_list[[j]]) == gene)
      if(length(idx) == 0) return(1)
      gene_de_list$de_list[[j]][idx, "p_val"]
    })
    stats::quantile(vec, probs = 0.75)
  })
  
  max(-log10(celltype_vec))
})
names(logpval_vec) <- Seurat::VariableFeatures(mbrain)
logpval_vec <- pmin(logpval_vec, 300)

rsquare_vec <- tiltedCCA:::postprocess_alignment(input_obj = multiSVD_obj,
                                                 bool_use_denoised = T,
                                                 seurat_obj = mbrain,
                                                 input_assay = 1,
                                                 seurat_assay = "RNA",
                                                 seurat_slot = "data")
logpval_vec <- logpval_vec[names(rsquare_vec)]
all(names(logpval_vec) == names(rsquare_vec))
stats::median(rsquare_vec[which(logpval_vec >= 1)])

png("../../../out/figures/main/10x_mouseembryo_differential_gene.png",
    height = 3500, width = 2500, res = 500, units = "px")
par(mar = c(5,5,4,1))
tiltedCCA:::plot_alignment(rsquare_vec = rsquare_vec,
                           logpval_vec = logpval_vec,
                           main = "Mouse Embryo E18 (10x, RNA+ATAC)\nGene differentiability vs. alignment",
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
                           lwd_polygon_bold = 5,
                           mark_median_xthres = 1)
graphics.off()
