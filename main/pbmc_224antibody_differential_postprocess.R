rm(list=ls())

library(Seurat)
library(Signac)
library(tiltedCCA)

load("../../../out/main/citeseq_bm25_preprocessed.RData")
ab_vec <- bm[["ADT"]]@var.features

load("../../../out/main/citeseq_pbmc224_differential.RData")
load("../../../out/main/citeseq_pbmc224_tiltedcca.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

##################3

Seurat::DefaultAssay(pbmc) <- "SCT"
gene_names <- Seurat::VariableFeatures(pbmc)
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
  names(celltype_vec) <- gene_de_list$level_vec
  # round(celltype_vec,2)
  
  max(-log10(celltype_vec))
})
names(logpval_vec) <- Seurat::VariableFeatures(pbmc)
logpval_vec <- pmin(logpval_vec, 300)

rsquare_vec <- tiltedCCA:::postprocess_modality_alignment(input_obj = multiSVD_obj,
                                                 bool_use_denoised = T,
                                                 seurat_obj = pbmc,
                                                 input_assay = 1,
                                                 min_subsample_cell = 5000,
                                                 seurat_assay = "SCT",
                                                 seurat_celltype_variable = "celltype.l2",
                                                 seurat_slot = "data",
                                                 verbose = 2)
all(names(logpval_vec) == names(rsquare_vec))
stats::median(rsquare_vec[which(logpval_vec >= 10)])

png("../../../out/figures/main/citeseq_pbmc224_differential_gene.png",
    height = 3500, width = 2500, res = 500, units = "px")
par(mar = c(5,5,4,1))
tiltedCCA:::plot_alignment(rsquare_vec = rsquare_vec,
                           logpval_vec = logpval_vec,
                           main = "PBMC (CITE-Seq, RNA+224 ADT)\nGene differentiability vs. alignment",
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
gene_depth <- log10(Matrix::rowSums(pbmc[["SCT"]]@counts[names(rsquare_vec),])+1)
gene_color <- sapply(gene_depth, function(x){
  col_palette_enrichment[which.min(abs(gene_breaks - x))]
})
ord_idx <- order(gene_depth, decreasing = F)

png("../../../out/figures/main/citeseq_pbmc224_differential_gene-colored.png",
    height = 3500, width = 2500, res = 500, units = "px")
# par(mar = c(5,5,4,1), bg = NA)
par(mar = c(5,5,4,1))
tiltedCCA:::plot_alignment(rsquare_vec = rsquare_vec[ord_idx],
                           logpval_vec = logpval_vec[ord_idx],
                           main = "PBMC (CITE-Seq, RNA+224 ADT)\nGene differentiability vs. alignment",
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

#################3

Cell_cycle <- c(cc.genes$s.genes[which(cc.genes$s.genes %in% gene_names)],
                cc.genes$g2m.genes[which(cc.genes$g2m.genes %in% gene_names)])

png("../../../out/figures/main/citeseq_pbmc224_differential_gene_Cell_cycle.png",
    height = 3500, width = 2500, res = 500, units = "px")
par(mar = c(5,5,4,1))
tiltedCCA:::plot_alignment(rsquare_vec = rsquare_vec,
                           logpval_vec = logpval_vec,
                           main = "PBMC (CITE-Seq, RNA+224 ADT)\nGene differentiability vs. alignment",
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

#########################
#########################
#########################

Seurat::DefaultAssay(pbmc) <- "ADT"
adt_names <- rownames(pbmc[["ADT"]]@counts)
logpval_vec <- sapply(1:length(adt_names), function(k){
  if(k %% floor(length(adt_names)/10) == 0) cat('*')
  adt <- adt_names[k]
  
  # cycle through all the celltypes
  celltype_vec <- sapply(1:length(adt_de_list$level_vec), function(i){
    idx <- which(adt_de_list$combn_mat == i, arr.ind = T)[,2]
    vec <-  sapply(idx, function(j){
      idx <- which(rownames(adt_de_list$de_list[[j]]) == adt)
      if(length(idx) == 0) return(1)
      adt_de_list$de_list[[j]][idx, "p_val"]
    })
    stats::quantile(vec, probs = 0.75)
  })
  names(celltype_vec) <- adt_de_list$level_vec
  # round(celltype_vec,2)
  
  max(-log10(celltype_vec))
})
names(logpval_vec) <- adt_names
logpval_vec <- pmin(logpval_vec, 300)

rsquare_vec <- tiltedCCA:::postprocess_modality_alignment(input_obj = multiSVD_obj,
                                                          bool_use_denoised = T,
                                                          seurat_obj = pbmc,
                                                          input_assay = 2,
                                                          min_subsample_cell = 5000,
                                                          seurat_assay = "ADT",
                                                          seurat_celltype_variable = "celltype.l2",
                                                          seurat_slot = "data",
                                                          verbose = 2)
all(names(logpval_vec) == names(rsquare_vec))
stats::median(rsquare_vec[which(logpval_vec >= 10)])

png("../../../out/figures/main/citeseq_pbmc224_differential_protein.png",
    height = 3500, width = 2500, res = 500, units = "px")
par(mar = c(5,5,4,1))
tiltedCCA:::plot_alignment(rsquare_vec = rsquare_vec,
                           logpval_vec = logpval_vec,
                           main = "PBMC (CITE-Seq, RNA+224 ADT)\nProtein differentiability vs. alignment",
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

#########
ab_vec_notin <- ab_vec[which(!ab_vec %in% rownames(pbmc[["ADT"]]))]
ab_vec_notin 
# [1] "CD11a"       "CD127-IL7Ra" "CD197-CCR7"  "CD278-ICOS"  "CD3"        
# [6] "CD38"        "CD4"         "CD56"        "HLA.DR"
sort(rownames(pbmc[["ADT"]]))
ab_vec_other <- c("CD11a/CD18", "CD127", "CD278", "CD3-1", 
                  "CD38-1", "CD38-2", "CD4-2", "CD56-1", 
                  "HLA-DR")
ab_vec_new <- c(ab_vec_other, ab_vec[which(ab_vec %in% rownames(pbmc[["ADT"]]))])
all(ab_vec_new %in% rownames(pbmc[["ADT"]]))

# adding jitter for visual clarity
logpval_vec2 <- logpval_vec
ab_idx <- which(names(logpval_vec2) %in% ab_vec_new)
idx <- ab_idx[which(logpval_vec2[ab_idx] >= 290)]
set.seed(10)
logpval_vec2[idx] <- logpval_vec2[idx] - runif(length(idx), min = 50, max = 200)

png("../../../out/figures/main/citeseq_pbmc224_differential_protein_highlight25.png",
    height = 3500, width = 2500, res = 500, units = "px")
par(mar = c(5,5,4,1))
tiltedCCA:::plot_alignment(rsquare_vec = rsquare_vec,
                           logpval_vec = logpval_vec2,
                           main = "PBMC (CITE-Seq, RNA+224 ADT)\nProtein differentiability vs. alignment",
                           bool_mark_ymedian = F,
                           bool_polygon_mean = T,
                           col_points = rgb(0.5, 0.5, 0.5, 0.15),
                           col_gene_highlight = rgb(82, 185, 44, maxColorValue = 255),
                           cex_axis = 1.5, 
                           cex_lab = 1.5,
                           cex_points = 2.5,
                           density = 10,
                           gene_names = ab_vec_new,
                           lty_polygon = 1,
                           lwd_grid_major = 2,
                           lwd_grid_minor = 1,
                           lwd_axis = 1.5,
                           lwd_axis_ticks = 1.5,
                           lwd_polygon = 2,
                           lwd_polygon_bold = 4)
graphics.off()

