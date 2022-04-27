rm(list=ls())
load("../../../out/main/citeseq_pbmc224_tiltedcca.RData")

library(Seurat); library(Signac)

input_obj <- multiSVD_obj
seurat_assay <- "SCT"

input_obj$param$svd_center_1 <- F
input_obj$param$svd_center_2 <- T
input_obj$param$svd_scale_1 <- F
input_obj$param$svd_scale_2 <- F

input_obj <- tiltedCCA:::.set_defaultAssay(input_obj, assay = 2)
dimred <- tiltedCCA:::.get_tCCAobj(input_obj, apply_postDimred = T, what = "distinct_mat")
seurat_umap <- Seurat::RunUMAP(dimred, 
                               assay = seurat_assay)
rownames(seurat_umap@cell.embeddings) <- rownames(dimred)
pbmc[["tmp"]] <- Seurat::CreateDimReducObject(seurat_umap@cell.embeddings,
                                              assay = seurat_assay)

plot1 <- Seurat::DimPlot(pbmc, reduction = "tmp",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster=FALSE)
plot1 <- plot1 + ggplot2::ggtitle(paste0("PBMC (CITE-Seq, RNA+224 ADT)\nCommon subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/citeseq_pbmc224_tcca-umap_distinct_version2.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

##########################

multiSVD_obj <- tiltedCCA:::tiltedCCA_decomposition(input_obj = multiSVD_obj,
                                                    verbose = 1,
                                                    bool_modality_1_full = F,
                                                    bool_modality_2_full = F)
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 2)
dimred <- tiltedCCA:::.get_tCCAobj(multiSVD_obj, apply_postDimred = T, what = "distinct_dimred")

seurat_assay <- "SCT"
seurat_umap <- Seurat::RunUMAP(dimred, 
                               assay = seurat_assay)
rownames(seurat_umap@cell.embeddings) <- rownames(dimred)
pbmc[["tmp"]] <- Seurat::CreateDimReducObject(seurat_umap@cell.embeddings,
                                              assay = seurat_assay)

plot1 <- Seurat::DimPlot(pbmc, reduction = "tmp",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster = FALSE)
plot1 <- plot1 + ggplot2::ggtitle(paste0("PBMC (CITE-Seq, RNA+224 ADT)\nCommon subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/citeseq_pbmc224_tcca-umap_distinct_version3.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")


