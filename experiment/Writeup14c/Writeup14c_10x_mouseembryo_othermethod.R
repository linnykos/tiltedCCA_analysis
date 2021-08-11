rm(list=ls()); set.seed(10)

library(Seurat); library(Signac); library(multiomicCCA)

load("../../../../out/Writeup14/Writeup14_10x_mouseembryo_preprocess.RData")
date_of_run <- Sys.time(); session_info <- sessionInfo()
metadata <- mbrain@meta.data

mat_1 <- mbrain[["RNA"]]@counts
Seurat::DefaultAssay(mbrain) <- "ATAC"
mat_2 <- mbrain[["ATAC"]]@counts


cell_idx <- which(metadata$label_Savercat %in% c("Oligodendrocyte", "Radial glia", 
                                                 "Forebrain GABAergic", "Neuroblast", 
                                                 "Glioblast", "Cortical or hippocampal glutamatergic"))

mbrain2 <- Seurat::CreateSeuratObject(counts = mat_1[,cell_idx])
mbrain2[["ATAC"]] <- Seurat::CreateAssayObject(counts = mat_2[,cell_idx])
mbrain2[["celltype"]] <- metadata$label_Savercat[cell_idx]

#############

set.seed(10)
Seurat::DefaultAssay(mbrain2) <- "RNA"
mbrain2 <- Seurat::SCTransform(mbrain2, verbose = T)
mbrain2 <- Seurat::RunPCA(mbrain2, verbose = F) 
mbrain2 <- Seurat::RunUMAP(mbrain2, dims = 1:50, reduction.name = 'umap.rna', 
                           reduction.key = 'rnaUMAP_')

set.seed(10)
Seurat::DefaultAssay(mbrain2) <- "ATAC"
mbrain2 <- Signac::RunTFIDF(mbrain2)
mbrain2 <- Signac::FindTopFeatures(mbrain2, min.cutoff = 'q0')
mbrain2 <- Signac::RunSVD(mbrain2)
mbrain2 <- Seurat::RunUMAP(mbrain2, reduction = 'lsi', dims = 2:50, 
                           reduction.name = "umap.atac", reduction.key = "atacUMAP_")

set.seed(10)
mbrain2 <- Seurat::FindMultiModalNeighbors(mbrain2, reduction.list = list("pca", "lsi"), 
                                           dims.list = list(1:50, 2:50))
mbrain2 <- Seurat::RunUMAP(mbrain2, nn.name = "weighted.nn", reduction.name = "wnn.umap", 
                           reduction.key = "wnnUMAP_")

plot1 <- Seurat::DimPlot(mbrain2, reduction = 'wnn.umap', group.by = 'celltype', 
                         label = TRUE, repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle("Mouse Embryo\nWNN") 
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14c/Writeup14c_10x_mouseembryo_wnn.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

#####################################

