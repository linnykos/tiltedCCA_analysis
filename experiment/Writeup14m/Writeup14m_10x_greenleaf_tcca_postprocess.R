rm(list=ls())
load("../../../../out/Writeup14m/Writeup14m_10x_greenleaf_geneactivity_tcca.RData")

library(Seurat)

plot1 <- Seurat::DimPlot(greenleaf, reduction = "common_tcca",
                         group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nCommon subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14m/Writeup14m_10x_greenleaf_geneactivity_umap_common-tcca.png"),
                plot1, device = "png", width = 5, height = 5, units = "in")

plot1 <-Seurat::DimPlot(greenleaf, reduction = "distinct1_tcca",
                        group.by = "celltype", label = TRUE,
                        repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nDistinct 1 subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14m/Writeup14m_10x_greenleaf_geneactivity_umap_distinct1-tcca.png"),
                plot1, device = "png", width = 5, height = 5, units = "in")

plot1 <-Seurat::DimPlot(greenleaf, reduction = "distinct2_tcca",
                        group.by = "celltype", label = TRUE,
                        repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nDistinct 2 subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14m/Writeup14m_10x_greenleaf_geneactivity_umap_distinct2-tcca.png"),
                plot1, device = "png", width = 5, height = 5, units = "in")

###############