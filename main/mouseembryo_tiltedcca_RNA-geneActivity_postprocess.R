rm(list=ls())
library(Seurat); library(Signac)

source("mouseembryo_colorPalette.R")
load("../../../out/main/10x_mouseembryo_tiltedcca_RNA-geneActivity.RData")

plot1 <- Seurat::DimPlot(mbrain, reduction = "common_tcca",
                         group.by = "label_Savercat", 
                         cols = col_palette)
plot1 <- plot1 + Seurat::NoLegend() + Seurat::NoAxes()
plot1 <- plot1 + ggplot2::ggtitle("")
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_mouseembryo_RNA-geneActivity_tcca-common-umap_cleaned.png"),
                plot1, device = "png", width = 3, height = 3, units = "in",
                dpi = 500)

plot1 <- Seurat::DimPlot(mbrain, reduction = "distinct1_tcca",
                         group.by = "label_Savercat", 
                         cols = col_palette)
plot1 <- plot1 + Seurat::NoLegend() + Seurat::NoAxes()
plot1 <- plot1 + ggplot2::ggtitle("")
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_mouseembryo_RNA-geneActivity_tcca-distinct1-umap_cleaned.png"),
                plot1, device = "png", width = 3, height = 3, units = "in",
                dpi = 500)

plot1 <- Seurat::DimPlot(mbrain, reduction = "distinct2_tcca",
                         group.by = "label_Savercat", 
                         cols = col_palette)
plot1 <- plot1 + Seurat::NoLegend() + Seurat::NoAxes()
plot1 <- plot1 + ggplot2::ggtitle("")
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_mouseembryo_RNA-geneActivity_tcca-distinct2-umap_cleaned.png"),
                plot1, device = "png", width = 3, height = 3, units = "in",
                dpi = 500)

############################

plot1 <- Seurat::DimPlot(mbrain, reduction = "common_tcca",
                         group.by = "label_Savercat", 
                         cols = col_palette)
plot1 <- plot1 + Seurat::NoLegend() + Seurat::NoAxes()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_mouseembryo_RNA-geneActivity_tcca-common-umap_cleaned.png"),
                plot1, device = "png", width = 3, height = 3, units = "in",
                dpi = 500)

set.seed(10)
mbrain <- Seurat::FindMultiModalNeighbors(mbrain, 
                                          reduction.list = list("pca", "pca.ga"), 
                                          weighted.nn.name = "weighted.nn2",
                                          dims.list = list(1:50, 2:50))
set.seed(10)
mbrain <- Seurat::RunUMAP(mbrain, nn.name = "weighted.nn2", reduction.name = "umap.wnn2", 
                          reduction.key = "wnnUMAP2_")

plot1 <- Seurat::DimPlot(mbrain, reduction = "umap.wnn2",
                         group.by = "label_Savercat", 
                         cols = col_palette)
plot1 <- plot1 + Seurat::NoLegend() + Seurat::NoAxes()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_mouseembryo_wnn-umap_RNA-geneActivity_cleaned.png"),
                plot1, device = "png", width = 3, height = 3, units = "in",
                dpi = 500)

# consensus pca
Seurat::DefaultAssay(mbrain) <- "SCT"
mat_1 <- Matrix::t(mbrain[["SCT"]]@data[Seurat::VariableFeatures(object = mbrain),])
Seurat::DefaultAssay(mbrain) <- "geneActivity"
mat_2 <- Matrix::t(mbrain[["geneActivity"]]@data[Seurat::VariableFeatures(object = mbrain),])

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
                                           dims_1 = 1:50, dims_2 = 2:50,
                                           dims_consensus = 1:49,
                                           center_1 = T, center_2 = T,
                                           recenter_1 = F, recenter_2 = F,
                                           rescale_1 = F, rescale_2 = F,
                                           scale_1 = T, scale_2 = T,
                                           verbose = 1)
consensus_dimred <- consensus_pca$dimred_consensus
colnames(consensus_dimred) <- paste0("consensusPCA2_", 1:ncol(consensus_dimred))
mbrain[["consensusPCA2"]] <- Seurat::CreateDimReducObject(consensus_dimred, 
                                                          assay = "SCT")

set.seed(10)
umap_res <- Seurat::RunUMAP(consensus_pca$dimred_consensus)
umap_mat <- umap_res@cell.embeddings
rownames(umap_mat) <- colnames(mbrain)
colnames(umap_mat) <- paste0("consensusUMAP2_", 1:ncol(umap_mat))
mbrain[["consensusUMAP2"]] <- Seurat::CreateDimReducObject(umap_mat, 
                                                           assay = "SCT")

plot1 <- Seurat::DimPlot(mbrain, reduction = "consensusUMAP2",
                         group.by = "label_Savercat", 
                         cols = col_palette)
plot1 <- plot1 + Seurat::NoLegend() + Seurat::NoAxes()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_mouseembryo_consensusPCA-umap_RNA-geneActivity_cleaned.png"),
                plot1, device = "png", width = 3, height = 3, units = "in",
                dpi = 500)


