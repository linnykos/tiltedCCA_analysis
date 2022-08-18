rm(list=ls())
library(Seurat)
library(MOFA2)

load("../../../out/main/citeseq_bm25_preprocessed.RData")
source("bm_25antibody_colorPalette.R")

model <- MOFA2::load_model("../../../out/main/citeseq_bm25_MOFA_model.hdf5")
factor_dims <- 1:MOFA2::get_dimensions(model)[["K"]]
factors <- MOFA2::get_factors(model, factors = factor_dims, groups = "all")
factors <- factors[[1]]

set.seed(10)
umap_res <- Seurat::RunUMAP(factors)
umap_mat <- umap_res@cell.embeddings
rownames(umap_mat) <- colnames(bm)
colnames(umap_mat) <- paste0("mofaUMAP_", 1:ncol(umap_mat))
bm[["mofaUMAP"]] <- Seurat::CreateDimReducObject(umap_mat, 
                                                    assay = "RNA")

plot1 <- Seurat::DimPlot(bm, reduction = "mofaUMAP",
                         group.by = "celltype.l2", 
                         cols = col_palette)
plot1 <- plot1 + Seurat::NoLegend() + Seurat::NoAxes()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_bm25_MOFA_cleaned.png"),
                plot1, device = "png", width = 3, height = 3, units = "in",
                dpi = 500)
