rm(list=ls())
load("../../../out/main/citeseq_bm25_preprocessed.RData")
source("bm_25antibody_colorPalette.R")

library(Seurat)
library(tiltedCCA)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

Seurat::DefaultAssay(bm) <- "RNA"
bm <- Seurat::ScaleData(bm) 
bm <- Seurat::RunPCA(bm, verbose = F)
Seurat::DefaultAssay(bm) <- "ADT"
bm <- Seurat::ScaleData(bm)
bm <- Seurat::RunPCA(bm, reduction.name = 'apca', verbose = F)

pca_rna <- bm[["pca"]]@cell.embeddings[,1:30]
pca_adt <- bm[["apca"]]@cell.embeddings[,1:18]
n <- nrow(pca_rna)
sing_rna <- apply(pca_rna, 2, function(x){sqrt(sum(x^2))})
sing_adt <- apply(pca_adt, 2, function(x){sqrt(sum(x^2))})
pca_rna <- tiltedCCA:::.mult_mat_vec(pca_rna, sqrt(n)/sing_rna)
pca_adt <- tiltedCCA:::.mult_mat_vec(pca_adt, sqrt(n)/sing_adt)

consensus_mat <- cbind(pca_rna, pca_adt)
set.seed(10)
umap_res <- Seurat::RunUMAP(consensus_mat)

umap_mat <- umap_res@cell.embeddings
rownames(umap_mat) <- colnames(bm)
bm[["consensusPCA"]] <- Seurat::CreateDimReducObject(umap_mat)

plot1 <- Seurat::DimPlot(bm, reduction = "consensusPCA",
                        group.by = "celltype.l2", label = TRUE,
                        repel = TRUE, label.size = 2.5,
                        cols = col_palette)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human BM (CITE-Seq, RNA+ADT)\nConsensus PCA"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_bm25_consensusPCA.png"),
                plot1, device = "png", width = 5.5, height = 5, units = "in")

plot1 <- Seurat::DimPlot(bm, reduction = "consensusPCA",
                         group.by = "celltype.l2", 
                         cols = col_palette)
plot1 <- plot1 + Seurat::NoLegend() + Seurat::NoAxes()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_bm25_consensusPCA_cleaned.png"),
                plot1, device = "png", width = 3, height = 3, units = "in",
                dpi = 500)

