rm(list=ls())
load("../../../out/main/citeseq_bm25_tcca.RData")
source("bm_25antibody_colorPalette.R")

plot1 <-Seurat::DimPlot(bm, reduction = "common_tcca",
                        group.by = "celltype.l2", label = TRUE,
                        repel = TRUE, label.size = 2.5,
                        cols = col_palette)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human BM (CITE-Seq, RNA+ADT)\nCommon subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_bm25_tcca-umap_common.png"),
                plot1, device = "png", width = 5.5, height = 5, units = "in")

plot1 <-Seurat::DimPlot(bm, reduction = "distinct1_tcca",
                        group.by = "celltype.l2", label = TRUE,
                        repel = TRUE, label.size = 2.5,
                        cols = col_palette)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human BM (CITE-Seq, RNA+ADT)\nRNA Distinct subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_bm25_tcca-umap_distinct1.png"),
                plot1, device = "png", width = 5.5, height = 5, units = "in")

plot1 <-Seurat::DimPlot(bm, reduction = "distinct2_tcca",
                        group.by = "celltype.l2", label = TRUE,
                        repel = TRUE, label.size = 2.5,
                        cols = col_palette)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human BM (CITE-Seq, RNA+ADT)\nADT Distinct subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_bm25_tcca-umap_distinct2.png"),
                plot1, device = "png", width = 5.5, height = 5, units = "in")

##########################

plot1 <- Seurat::DimPlot(bm, reduction = "rna.umap",
                         group.by = "celltype.l2", 
                         cols = col_palette)
plot1 <- plot1 + Seurat::NoLegend() + Seurat::NoAxes()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_bm25_rna-umap_cleaned.png"),
                plot1, device = "png", width = 3, height = 3, units = "in",
                dpi = 500)

plot1 <- Seurat::DimPlot(bm, reduction = "adt.umap",
                         group.by = "celltype.l2", 
                         cols = col_palette)
plot1 <- plot1 + Seurat::NoLegend() + Seurat::NoAxes()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_bm25_adt-umap_cleaned.png"),
                plot1, device = "png", width = 3, height = 3, units = "in",
                dpi = 500)

plot1 <- Seurat::DimPlot(bm, reduction = "common_tcca",
                         group.by = "celltype.l2", 
                         cols = col_palette)
plot1 <- plot1 + Seurat::NoLegend() + Seurat::NoAxes()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_bm25_tcca-umap_common_cleaned.png"),
                plot1, device = "png", width = 3, height = 3, units = "in",
                dpi = 500)

plot1 <- Seurat::DimPlot(bm, reduction = "distinct1_tcca",
                         group.by = "celltype.l2", 
                         cols = col_palette)
plot1 <- plot1 + Seurat::NoLegend() + Seurat::NoAxes()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_bm25_tcca-umap_distinct1_cleaned.png"),
                plot1, device = "png", width = 3, height = 3, units = "in",
                dpi = 500)

plot1 <- Seurat::DimPlot(bm, reduction = "distinct2_tcca",
                         group.by = "celltype.l2", 
                         cols = col_palette)
plot1 <- plot1 + Seurat::NoLegend() + Seurat::NoAxes()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_bm25_tcca-umap_distinct2_cleaned.png"),
                plot1, device = "png", width = 3, height = 3, units = "in",
                dpi = 500)