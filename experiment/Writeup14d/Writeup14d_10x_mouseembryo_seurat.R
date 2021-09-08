rm(list=ls()); set.seed(10)

library(Seurat); library(Signac)

load("../../../../out/Writeup14/Writeup14_10x_mouseembryo_preprocess.RData")
date_of_run <- Sys.time(); session_info <- sessionInfo()

metadata <- mbrain@meta.data
cell_idx <- which(metadata$label_Savercat %in% c("Oligodendrocyte", "Radial glia", 
                                                 "Forebrain GABAergic", "Neuroblast", 
                                                 "Glioblast", "Cortical or hippocampal glutamatergic"))

n <- nrow(metadata)
keep_vec <- rep(0, n)
keep_vec[cell_idx] <- 1
mbrain[["keep"]] <- keep_vec

mbrain2 <- subset(mbrain, keep == 1)

##########

mbrain2[["umap"]] <- NULL
mbrain2[["umap.atac"]] <- NULL
mbrain2[["umap.wnn"]] <- NULL

Seurat::DefaultAssay(mbrain) <- "SCT"
set.seed(10)
mbrain2 <- Seurat::RunUMAP(mbrain2, reduction = 'pca', dims = 1:30, assay = 'SCT', 
              reduction.name = 'umap', reduction.key = 'sctUMAP_')
set.seed(10)
mbrain2 <- RunUMAP(mbrain2, reduction = 'lsi', dims = 2:50, assay = 'ATAC', 
              reduction.name = 'umap.atac', reduction.key = 'atacUMAP_')

mbrain2[["SCT.weight"]] <- NULL
mbrain2[["ATAC.weight"]] <- NULL
set.seed(10)
mbrain2 <- Seurat::FindMultiModalNeighbors(
  mbrain2,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:30, 2:50))
set.seed(10)
mbrain2 <- Seurat::RunUMAP(mbrain2, 
                           nn.name = "weighted.nn", 
                           reduction.name = "wnn.umap", 
                           reduction.key = "wnnUMAP_")
set.seed(10)
mbrain2 <- Seurat::FindClusters(mbrain2, 
                                graph.name = "wsnn", 
                                algorithm = 3, 
                                resolution = 2, 
                                verbose = T)

save(mbrain2, file = "../../../../out/Writeup14d/Writeup14d_10x_mouseembryo_seurat.RData")

######################

plot1 <- Seurat::DimPlot(mbrain2, reduction = "umap", 
                         group.by = "label_Savercat", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle("Mouse Embryo (RNA)\nSeurat baseline")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = "../../../../out/figures/Writeup14d/Writeup14d_10x_mouseembryo_seurat_rna_umap2.png",
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <- Seurat::DimPlot(mbrain2, reduction = "umap.atac", 
                         group.by = "label_Savercat", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle("Mouse Embryo (ATAC)\nSeurat baseline")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = "../../../../out/figures/Writeup14d/Writeup14d_10x_mouseembryo_seurat_atac_umap2.png",
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <- Seurat::DimPlot(mbrain2, reduction = "wnn.umap", 
                         group.by = "label_Savercat", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle("Mouse Embryo (WNN)\nSeurat baseline")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = "../../../../out/figures/Writeup14d/Writeup14d_10x_mouseembryo_seurat_wnn_umap2.png",
                plot1, device = "png", width = 6, height = 5, units = "in")

###

plot1 <- Seurat::DimPlot(mbrain2, reduction = "umap", 
                         group.by = "wsnn_res.2", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle("Mouse Embryo (RNA)\nSeurat baseline")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = "../../../../out/figures/Writeup14d/Writeup14d_10x_mouseembryo_seurat_rna_umap2_clustering.png",
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <- Seurat::DimPlot(mbrain2, reduction = "umap.atac", 
                         group.by = "wsnn_res.2", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle("Mouse Embryo (ATAC)\nSeurat baseline")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = "../../../../out/figures/Writeup14d/Writeup14d_10x_mouseembryo_seurat_atac_umap2_clustering.png",
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <- Seurat::DimPlot(mbrain2, reduction = "wnn.umap", 
                         group.by = "wsnn_res.2", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle("Mouse Embryo (WNN)\nSeurat baseline")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = "../../../../out/figures/Writeup14d/Writeup14d_10x_mouseembryo_seurat_wnn_umap2_clustering.png",
                plot1, device = "png", width = 6, height = 5, units = "in")


