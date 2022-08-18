rm(list=ls())
library(Seurat)

load("../../../out/main/citeseq_bm25_JIVE.RData")
source("bm_25antibody_colorPalette.R")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

####################

set.seed(10)
umap_res <- Seurat::RunUMAP(scai_res$H)
umap_mat <- umap_res@cell.embeddings
rownames(umap_mat) <- colnames(bm)
colnames(umap_mat) <- paste0("scaiUMAP_", 1:ncol(umap_mat))
bm[["scaiUMAP"]] <- Seurat::CreateDimReducObject(umap_mat, 
                                                 assay = "RNA")

plot1 <- Seurat::DimPlot(bm, reduction = "scaiUMAP",
                         group.by = "celltype.l2", 
                         cols = col_palette)
plot1 <- plot1 + Seurat::NoLegend() + Seurat::NoAxes()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_bm25_scAI_cleaned.png"),
                plot1, device = "png", width = 3, height = 3, units = "in",
                dpi = 500)
