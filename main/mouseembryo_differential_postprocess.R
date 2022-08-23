rm(list=ls())
load("../../../out/main/10x_mouseembryo_differential.RData")
load("../../../out/main/10x_mouseembryo_tiltedcca.RData")

library(Seurat)
library(Signac)
library(tiltedCCA)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

#########

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

rsquare_vec <- tiltedCCA:::postprocess_modality_alignment(input_obj = multiSVD_obj,
                                                          bool_use_denoised = T,
                                                          seurat_obj = mbrain,
                                                          input_assay = 1,
                                                          seurat_assay = "RNA",
                                                          seurat_slot = "data")
logpval_vec <- logpval_vec[names(rsquare_vec)]
all(names(logpval_vec) == names(rsquare_vec))
stats::median(rsquare_vec[which(logpval_vec >= 10)])

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
                           mark_median_xthres = 10)
graphics.off()

col_palette_enrichment <- paste0(grDevices::colorRampPalette(c('lightgrey', 'blue'))(100), "33")
gene_breaks <- seq(1, 5.5, length.out = 100)
gene_depth <- log10(Matrix::rowSums(mbrain[["RNA"]]@counts[names(rsquare_vec),])+1)
gene_color <- sapply(gene_depth, function(x){
  col_palette_enrichment[which.min(abs(gene_breaks - x))]
})
ord_idx <- order(gene_depth, decreasing = F)

png("../../../out/figures/main/10x_mouseembryo_differential_gene-colored.png",
    height = 3500, width = 2500, res = 500, units = "px")
# par(mar = c(5,5,4,1), bg = NA)
par(mar = c(5,5,4,1))
tiltedCCA:::plot_alignment(rsquare_vec = rsquare_vec[ord_idx],
                           logpval_vec = logpval_vec[ord_idx],
                           main = "Mouse Embryo E18 (10x, RNA+ATAC)\nGene differentiability vs. alignment",
                           bool_mark_ymedian = T,
                           col_gene_highlight_border = rgb(255, 205, 87, 255*0.5, maxColorValue = 255),
                           col_points = gene_color[ord_idx],
                           cex_axis = 1.5, 
                           cex_lab = 1.5,
                           cex_points = 2.5,
                           lty_polygon = 2,
                           lwd_grid_major = 2,
                           lwd_grid_minor = 1,
                           lwd_axis = 1.5,
                           lwd_axis_ticks = 1.5,
                           lwd_polygon_bold = 5,
                           mark_median_xthres = 10)
graphics.off()

###################################

Cell_cycle <- c(cc.genes$s.genes[which(cc.genes$s.genes %in% gene_names)],
                cc.genes$g2m.genes[which(cc.genes$g2m.genes %in% gene_names)])

png(paste0("../../../out/figures/main/10x_mouseembryo_differential_gene_Cell_cycle.png"),
    height = 3500, width = 2500, res = 500, units = "px")
par(mar = c(5,5,4,1))
tiltedCCA:::plot_alignment(rsquare_vec = rsquare_vec,
                           logpval_vec = logpval_vec,
                           main = "Mouse Embryo E18 (10x, RNA+ATAC)\nGene differentiability vs. alignment",
                           bool_mark_ymedian = F,
                           bool_polygon_mean = T,
                           col_points = rgb(0.5, 0.5, 0.5, 0.1),
                           col_gene_highlight = "black",
                           cex_axis = 1.5, 
                           cex_lab = 1.5,
                           cex_points = 2.5,
                           density = 10,
                           gene_names = Cell_cycle,
                           lty_polygon = 1,
                           lwd_grid_major = 2,
                           lwd_grid_minor = 1,
                           lwd_axis = 1.5,
                           lwd_axis_ticks = 1.5,
                           lwd_polygon = 2,
                           lwd_polygon_bold = 4)
graphics.off()
