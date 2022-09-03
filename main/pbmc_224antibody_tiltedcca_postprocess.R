rm(list=ls())
library(Seurat); library(Signac)

load("../../../out/main/citeseq_pbmc224_tiltedcca.RData")
source("pbmc_224antibody_colorPalette.R")

plot1 <- Seurat::DimPlot(pbmc, reduction = "common_tcca",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         cols = col_palette,
                         raster = FALSE)
plot1 <- plot1 + ggplot2::ggtitle(paste0("PBMC (CITE-Seq, RNA+224 ADT)\nCommon subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_pbmc224_tcca-umap_common.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot2 <- Seurat::DimPlot(pbmc, reduction = "distinct1_tcca",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         cols = col_palette,
                         raster = FALSE)
plot2 <- plot2 + ggplot2::ggtitle(paste0("PBMC (CITE-Seq, RNA+224 ADT)\nRNA distinct subspace"))
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_pbmc224_tcca-umap_distinct1.png"),
                plot2, device = "png", width = 6, height = 5, units = "in")

plot3 <- Seurat::DimPlot(pbmc, reduction = "distinct2_tcca",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         cols = col_palette,
                         raster = FALSE)
plot3 <- plot3 + ggplot2::ggtitle(paste0("PBMC (CITE-Seq, RNA+224 ADT)\nADT distinct subspace"))
plot3 <- plot3 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_pbmc224_tcca-umap_distinct2.png"),
                plot3, device = "png", width = 6, height = 5, units = "in")

###############################

plot1 <- Seurat::DimPlot(pbmc, reduction = "common_tcca",
                         group.by = "celltype.l2",
                         cols = col_palette,
                         raster = FALSE)
plot1 <- plot1 + Seurat::NoLegend() + Seurat::NoAxes()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_pbmc224_tcca-umap_common_cleaned.png"),
                plot1, device = "png", width = 3, height = 3, units = "in",
                dpi = 500)

plot2 <- Seurat::DimPlot(pbmc, reduction = "distinct1_tcca",
                         group.by = "celltype.l2", 
                         cols = col_palette,
                         raster = FALSE)
plot2 <- plot2 + Seurat::NoLegend() + Seurat::NoAxes()
plot2 <- plot2 + ggplot2::ggtitle("")
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_pbmc224_tcca-umap_distinct1_cleaned.png"),
                plot2, device = "png", width = 3, height = 3, units = "in",
                dpi = 500)

plot3 <- Seurat::DimPlot(pbmc, reduction = "distinct2_tcca",
                         group.by = "celltype.l2",
                         cols = col_palette,
                         raster = FALSE)
plot3 <- plot3 + Seurat::NoLegend() + Seurat::NoAxes()
plot3 <- plot3 + ggplot2::ggtitle("")
plot3 <- plot3 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_pbmc224_tcca-umap_distinct2_cleaned.png"),
                plot3, device = "png", width = 3, height = 3, units = "in",
                dpi = 500)

plot3 <- Seurat::DimPlot(pbmc, reduction = "wnn.umap",
                         group.by = "celltype.l2",
                         cols = col_palette,
                         raster = FALSE)
plot3 <- plot3 + Seurat::NoLegend() + Seurat::NoAxes()
plot3 <- plot3 + ggplot2::ggtitle("")
plot3 <- plot3 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_pbmc224_wnn-umap_cleaned.png"),
                plot3, device = "png", width = 3, height = 3, units = "in",
                dpi = 500)

###########################

# consensus pca
Seurat::DefaultAssay(pbmc) <- "SCT"
mat_1 <- Matrix::t(pbmc[["SCT"]]@data[Seurat::VariableFeatures(object = pbmc),])
Seurat::DefaultAssay(pbmc) <- "ADT"
mat_2 <- Matrix::t(pbmc[["ADT"]]@data)

mat_1b <- mat_1
sd_vec <- sparseMatrixStats::colSds(mat_1b)
if(any(sd_vec <= 1e-6)){
  mat_1b <- mat_1b[,-which(sd_vec <= 1e-6)]
}

mat_2b <- mat_2
sd_vec <- sparseMatrixStats::colSds(mat_2b)
if(any(sd_vec <= 1e-6)){
  mat_2b <- mat_2b[,-which(sd_vec <= 1e-6)]
}

consensus_pca <- tiltedCCA:::consensus_pca(mat_1 = mat_1b, mat_2 = mat_2b,
                                           dims_1 = 1:40, dims_2 = 1:25,
                                           dims_consensus = 1:40,
                                           center_1 = T, center_2 = T,
                                           recenter_1 = F, recenter_2 = F,
                                           rescale_1 = F, rescale_2 = F,
                                           scale_1 = T, scale_2 = T,
                                           verbose = 1)
consensus_dimred <- consensus_pca$dimred_consensus
colnames(consensus_dimred) <- paste0("consensusPCA_", 1:ncol(consensus_dimred))
pbmc[["consensusPCA"]] <- Seurat::CreateDimReducObject(consensus_dimred, 
                                                       assay = "SCT")

set.seed(10)
umap_res <- Seurat::RunUMAP(consensus_pca$dimred_consensus)
umap_mat <- umap_res@cell.embeddings
rownames(umap_mat) <- colnames(pbmc)
colnames(umap_mat) <- paste0("consensusUMAP_", 1:ncol(umap_mat))
pbmc[["consensusUMAP"]] <- Seurat::CreateDimReducObject(umap_mat, 
                                                        assay = "SCT")

plot3 <- Seurat::DimPlot(pbmc, reduction = "consensusUMAP",
                         group.by = "celltype.l2",
                         cols = col_palette,
                         raster = FALSE)
plot3 <- plot3 + Seurat::NoLegend() + Seurat::NoAxes()
plot3 <- plot3 + ggplot2::ggtitle("")
plot3 <- plot3 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_pbmc224_consensuPCA-umap_cleaned.png"),
                plot3, device = "png", width = 3, height = 3, units = "in",
                dpi = 500)

