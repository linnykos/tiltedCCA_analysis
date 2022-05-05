rm(list=ls())
load("../../../out/main/citeseq_pbmc224_tiltedcca-2000genes.RData")
source("pbmc_224antibody_colorPalette.R")

library(Seurat); library(Signac)

plot1 <- Seurat::DimPlot(pbmc, reduction = "common_tcca",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         cols = col_palette,
                         raster = FALSE)
plot1 <- plot1 + ggplot2::ggtitle(paste0("PBMC (CITE-Seq, RNA+224 ADT)\nCommon subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_pbmc224-2000genes_tcca-umap_common.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot2 <- Seurat::DimPlot(pbmc, reduction = "distinct1_tcca",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         cols = col_palette,
                         raster = FALSE)
plot2 <- plot2 + ggplot2::ggtitle(paste0("PBMC (CITE-Seq, RNA+224 ADT)\nRNA distinct subspace"))
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_pbmc224-2000genes_tcca-umap_distinct1.png"),
                plot2, device = "png", width = 6, height = 5, units = "in")

plot3 <- Seurat::DimPlot(pbmc, reduction = "distinct2_tcca",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         cols = col_palette,
                         raster = FALSE)
plot3 <- plot3 + ggplot2::ggtitle(paste0("PBMC (CITE-Seq, RNA+224 ADT)\nADT distinct subspace"))
plot3 <- plot3 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_pbmc224-2000genes_tcca-umap_distinct2.png"),
                plot3, device = "png", width = 6, height = 5, units = "in")

###############################

load("../../../out/main/citeseq_pbmc224_differential.RData")
load("../../../out/main/citeseq_pbmc224_tiltedcca-2000genes.RData")

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

rsquare_vec <- tiltedCCA:::postprocess_alignment(input_obj = multiSVD_obj,
                                                 bool_use_denoised = T,
                                                 seurat_obj = pbmc,
                                                 input_assay = 1,
                                                 min_subsample_cell = 5000,
                                                 seurat_assay = "SCT",
                                                 seurat_celltype_variable = "celltype.l2",
                                                 seurat_slot = "data",
                                                 verbose = 2)
logpval_vec <- logpval_vec[names(rsquare_vec)]
all(names(logpval_vec) == names(rsquare_vec))
stats::median(rsquare_vec[which(logpval_vec >=1)])
quantile(logpval_vec)
max_val <- max(logpval_vec[which(!is.infinite(logpval_vec))])
logpval_vec <- pmin(logpval_vec, 300)

png("../../../out/figures/main/citeseq_pbmc224-2000genes_differential_gene.png",
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
                           lwd_polygon_bold = 5)
graphics.off()

