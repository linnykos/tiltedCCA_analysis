rm(list=ls())
library(Seurat)
library(Signac)

load("../../../../../Public/10x_trevino_simplified.RData")
seurat_obj

table(seurat_obj$celltype)

set.seed(10)
Seurat::DefaultAssay(seurat_obj) <- "SCT"
seurat_obj <- Seurat::ScaleData(seurat_obj)
seurat_obj <- Seurat::RunPCA(seurat_obj, 
                             verbose = FALSE)
seurat_obj <- Seurat::RunUMAP(seurat_obj, 
                              dims = 1:50, 
                              reduction.name="umap.rna", 
                              reduction.key="rnaUMAP_")

set.seed(10)
seurat_obj <-  Signac::RunSVD(seurat_obj)  
seurat_obj <- Seurat::RunUMAP(seurat_obj,
                              reduction="lsi", 
                              dims=2:50, 
                              reduction.name="umap.atac", 
                              reduction.key="atacUMAP_")

seurat_obj <- tiltedCCA::rotate_seurat_embeddings(
  seurat_obj = seurat_obj,
  source_embedding = "umap.rna",
  target_embedding = "umap.atac"
)

####

col_palette <- c(
  "Cyc. Prog." = rgb(213, 163, 98, maxColorValue = 255),
  "GluN2" = rgb(122, 179, 232, maxColorValue = 255),
  "GluN3"= rgb(174, 198, 235, maxColorValue = 255),
  "GluN4" = rgb(217, 227, 132, maxColorValue = 255),
  "GluN5" = rgb(127, 175, 123, maxColorValue = 255),
  "nIPC/GluN1" = rgb(114, 169, 158, maxColorValue = 255),
  "RG" = rgb(197, 125, 95, maxColorValue = 255)
)


p1 <- Seurat::DimPlot(seurat_obj, 
                      reduction = "umap.rna",
                      group.by = "celltype", 
                      cols = col_palette, 
                      label = T, repel = T, label.size = 2.5)
p1 <- p1 + Seurat::NoLegend()
p1 <- p1 + ggplot2::ggtitle("RNA") + ggplot2::labs(x = "", y = "")
p1 <- p1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

p2 <- Seurat::DimPlot(seurat_obj, 
                      reduction = "umap.atac",
                      group.by = "celltype", 
                      cols = col_palette, 
                      label = T, repel = T, label.size = 2.5)
p2 <- p2 + Seurat::NoLegend()
p2 <- p2 + ggplot2::ggtitle("ATAC") + ggplot2::labs(x = "", y = "")
p2 <- p2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

p_all <- cowplot::plot_grid(p1, p2, ncol = 2)
ggplot2::ggsave(filename = "vignettes/embryo_umap.png",
                p_all, device = "png", height = 1250, width = 2250, units = "px",
                dpi = 300)

#######################

Seurat::DefaultAssay(seurat_obj) <- "SCT"
mat_1 <- Matrix::t(seurat_obj[["SCT"]]@data)
Seurat::DefaultAssay(seurat_obj) <- "ATAC"
mat_2 <- Matrix::t(seurat_obj[["ATAC"]]@data)

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

set.seed(10)
multiSVD_obj <- tiltedCCA:::create_multiSVD(mat_1 = mat_1b, mat_2 = mat_2b,
                                            dims_1 = 1:50, dims_2 = 2:50,
                                            center_1 = T, center_2 = F,
                                            normalize_row = T,
                                            normalize_singular_value = T,
                                            recenter_1 = F, recenter_2 = T,
                                            rescale_1 = F, rescale_2 = T,
                                            scale_1 = T, scale_2 = F,
                                            verbose = 1)
multiSVD_obj <- tiltedCCA:::form_metacells(input_obj = multiSVD_obj,
                                           large_clustering_1 = NULL, 
                                           large_clustering_2 = NULL, 
                                           num_metacells = NULL,
                                           verbose = 1)
multiSVD_obj <- tiltedCCA:::compute_snns(input_obj = multiSVD_obj,
                                         latent_k = 20,
                                         num_neigh = 15,
                                         bool_cosine = T,
                                         bool_intersect = F,
                                         min_deg = 15,
                                         verbose = 1)
multiSVD_obj2 <- tiltedCCA:::tiltedCCA(input_obj = multiSVD_obj,
                                       verbose = 1)
multiSVD_obj2 <- tiltedCCA:::fine_tuning(input_obj = multiSVD_obj2,
                                         verbose = 1)
multiSVD_obj2 <- tiltedCCA:::tiltedCCA_decomposition(input_obj = multiSVD_obj2,
                                                     verbose = 1,
                                                     bool_modality_1_full = T,
                                                     bool_modality_2_full = F)

set.seed(10)
seurat_obj[["common_tcca"]] <- tiltedCCA::create_SeuratDim(input_obj = multiSVD_obj2,
                                                           what = "common",
                                                           aligned_umap_assay = "umap.rna",
                                                           seurat_obj = seurat_obj,
                                                           seurat_assay = "SCT",
                                                           verbose = 1)
set.seed(10)
seurat_obj[["distinct1_tcca"]] <- tiltedCCA::create_SeuratDim(input_obj = multiSVD_obj2,
                                                              what = "distinct_1",
                                                              aligned_umap_assay = "umap.rna",
                                                              seurat_obj = seurat_obj,
                                                              seurat_assay = "SCT",
                                                              verbose = 1)
set.seed(10)
seurat_obj[["distinct2_tcca"]] <- tiltedCCA::create_SeuratDim(input_obj = multiSVD_obj2,
                                                              what = "distinct_2",
                                                              aligned_umap_assay = "umap.rna",
                                                              seurat_obj = seurat_obj,
                                                              seurat_assay = "SCT",
                                                              verbose = 1)

p1 <- Seurat::DimPlot(seurat_obj, 
                      reduction = "common_tcca",
                      group.by = "celltype", 
                      cols = col_palette, 
                      label = T, repel = T, label.size = 2.5)
p1 <- p1 + Seurat::NoLegend()
p1 <- p1 + ggplot2::ggtitle("Common") + ggplot2::labs(x = "", y = "")
p1 <- p1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

p2 <- Seurat::DimPlot(seurat_obj, 
                      reduction = "distinct1_tcca",
                      group.by = "celltype", 
                      cols = col_palette, 
                      label = T, repel = T, label.size = 2.5)
p2 <- p2 + Seurat::NoLegend()
p2 <- p2 + ggplot2::ggtitle("RNA Distinct") + ggplot2::labs(x = "", y = "")
p2 <- p2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

p3 <- Seurat::DimPlot(seurat_obj, 
                      reduction = "distinct2_tcca",
                      group.by = "celltype", 
                      cols = col_palette, 
                      label = T, repel = T, label.size = 2.5)
p3 <- p3 + Seurat::NoLegend()
p3 <- p3 + ggplot2::ggtitle("ATAC Distinct") + ggplot2::labs(x = "", y = "")
p3 <- p3 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

p_all <- cowplot::plot_grid(p1, p2, p3, ncol = 3)
ggplot2::ggsave(filename = "vignettes/embryo_tiltedcca.png",
                p_all, device = "png", height = 1250, width = 3500, units = "px",
                dpi = 300)
