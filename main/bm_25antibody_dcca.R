rm(list=ls())
load("../../../out/main/citeseq_bm25_tcca.RData")
source("bm_25antibody_colorPalette.R")

library(Seurat); library(Signac)

multiSVD_obj2 <- tiltedCCA:::tiltedCCA(input_obj = multiSVD_obj,
                                       fix_tilt_perc = 0.5,
                                       verbose = 1)
multiSVD_obj2 <- tiltedCCA:::tiltedCCA_decomposition(input_obj = multiSVD_obj2,
                                                     verbose = 1)
set.seed(10)
bm[["common_tcca2"]] <- tiltedCCA:::create_SeuratDim(input_obj = multiSVD_obj2,
                                                     what = "common",
                                                     aligned_umap_assay = "rna.umap",
                                                     seurat_obj = bm,
                                                     seurat_assay = "RNA",
                                                     verbose = 1)
set.seed(10)
bm[["distinct1_tcca2"]] <- tiltedCCA:::create_SeuratDim(input_obj = multiSVD_obj2,
                                                     what = "distinct_1",
                                                     aligned_umap_assay = "rna.umap",
                                                     seurat_obj = bm,
                                                     seurat_assay = "RNA",
                                                     verbose = 1)

set.seed(10)
bm[["distinct2_tcca2"]] <- tiltedCCA:::create_SeuratDim(input_obj = multiSVD_obj2,
                                                        what = "distinct_2",
                                                        aligned_umap_assay = "rna.umap",
                                                        seurat_obj = bm,
                                                        seurat_assay = "RNA",
                                                        verbose = 1)

plot1 <- Seurat::DimPlot(bm, reduction = "common_tcca2",
                         group.by = "celltype.l2", 
                         cols = col_palette)
plot1 <- plot1 + Seurat::NoLegend() + Seurat::NoAxes()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_bm25_dcca-common_cleaned.png"),
                plot1, device = "png", width = 3, height = 3, units = "in",
                dpi = 500)

plot1 <- Seurat::DimPlot(bm, reduction = "distinct1_tcca2",
                         group.by = "celltype.l2", 
                         cols = col_palette)
plot1 <- plot1 + Seurat::NoLegend() + Seurat::NoAxes()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_bm25_dcca-distinct1_cleaned.png"),
                plot1, device = "png", width = 3, height = 3, units = "in",
                dpi = 500)

plot1 <- Seurat::DimPlot(bm, reduction = "distinct2_tcca2",
                         group.by = "celltype.l2", 
                         cols = col_palette)
plot1 <- plot1 + Seurat::NoLegend() + Seurat::NoAxes()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_bm25_dcca-distinct2_cleaned.png"),
                plot1, device = "png", width = 3, height = 3, units = "in",
                dpi = 500)