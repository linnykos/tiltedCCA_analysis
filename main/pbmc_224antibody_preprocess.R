rm(list=ls())
load("../../../data/CITE_PBMC_RNA-228Protein/pbmc.RData")
library(Seurat)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

#######################

pbmc[["spca"]] <- NULL
pbmc[["umap"]] <- NULL
pbmc[["aumap"]] <- NULL
pbmc[["wknn"]] <- NULL
pbmc[["wsnn"]] <- NULL
pbmc[["wnn.umap"]] <- NULL

pbmc[["SCT"]]@scale.data <- tcrossprod(pbmc[["pca"]]@feature.loadings, 
                                       pbmc[["pca"]]@cell.embeddings)
pbmc[["ADT"]]@scale.data <- tcrossprod(pbmc[["apca"]]@feature.loadings, 
                                       pbmc[["apca"]]@cell.embeddings)
dim(pbmc[["SCT"]]@scale.data)
dim(pbmc[["ADT"]]@scale.data)

set.seed(10)
pbmc <- Seurat::RunUMAP(pbmc, reduction = 'pca', dims = 1:40, assay = 'SCT',
                        reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')

set.seed(10)
pbmc <- Seurat::RunUMAP(pbmc, reduction = 'apca', dims = 1:50, assay = 'ADT',
                        reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')

# the following takes a hot minute with 160K cells
set.seed(10)
pbmc <- Seurat::FindMultiModalNeighbors(pbmc, 
                                        reduction.list = list("pca", "apca"),
                                        dims.list = list(1:40, 1:50), 
                                        modality.weight.name = "RNA.weight")
set.seed(10)
pbmc <- Seurat::RunUMAP(pbmc, 
                        nn.name = "weighted.nn", 
                        reduction.name = "wnn.umap", 
                        reduction.key = "wnnUMAP_")

save(pbmc, date_of_run, session_info,
     file = "../../../out/main/citeseq_pbmc224_preprocessed.RData")

###################

plot1 <- Seurat::DimPlot(pbmc, reduction = "rna.umap",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster=FALSE)
plot1 <- plot1 + ggplot2::ggtitle(paste0("PBMC (CITE-Seq, RNA+224 ADT)\nRNA UMAP"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_pbmc224_rna-umap.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot2 <- Seurat::DimPlot(pbmc, reduction = "adt.umap",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster=FALSE)
plot2 <- plot2 + ggplot2::ggtitle(paste0("PBMC (CITE-Seq, RNA+224 ADT)\nADT UMAP"))
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_pbmc224_adt-umap.png"),
                plot2, device = "png", width = 6, height = 5, units = "in")

plot3 <- Seurat::DimPlot(pbmc, reduction = "wnn.umap",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster=FALSE)
plot3 <- plot3 + ggplot2::ggtitle(paste0("PBMC (CITE-Seq, RNA+224 ADT)\nWNN UMAP"))
plot3 <- plot3 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_pbmc224_wnn-umap.png"),
                plot3, device = "png", width = 6, height = 5, units = "in")
