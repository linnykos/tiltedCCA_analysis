rm(list=ls())
source("bm_97antibody_colorPalette.R")

library(Seurat)
library(dplyr)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

bm <- readRDS("~/nzhanglab/data/Triana_nature_bonemarrow/WTA_projected.rds")
bm$ct <- factor(bm$ct, levels = sort(levels(bm$ct)))

plot1 <-Seurat::DimPlot(bm, reduction = "bothUMAP",
                        group.by = "ct", label = TRUE,
                        repel = TRUE, label.size = 2.5,
                        cols = col_palette)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human BM (Abseq, RNA+ADT)\nSupplied: bothUMAP"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/abseq_bm97_supplied-bothUMAP.png"),
                plot1, device = "png", width = 9, height = 5, units = "in")

plot1 <-Seurat::DimPlot(bm, reduction = "MOFAUMAP",
                        group.by = "ct", label = TRUE,
                        repel = TRUE, label.size = 2.5,
                        cols = col_palette)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human BM (Abseq, RNA+ADT)\nSupplied: MOFAUMAP"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/abseq_bm97_supplied-MOFAUMAP.png"),
                plot1, device = "png", width = 9, height = 5, units = "in")

##################

DefaultAssay(bm) <- "RNA"
bm <- Seurat::RunPCA(bm, verbose = F)

DefaultAssay(bm) <- "AB"
bm[["AB"]]@var.features <- rownames(bm)
bm <- Seurat::RunPCA(bm, reduction.name = 'apca',
                     verbose = F)

set.seed(10)
bm <- Seurat::RunUMAP(bm, reduction = 'pca', dims = 1:30, assay = 'RNA',
                      reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
set.seed(10)
bm <- Seurat::RunUMAP(bm, 
                      reduction = 'apca', 
                      dims = 1:30, assay = 'ADT',
                      reduction.name = 'adt.umap', 
                      reduction.key = 'adtUMAP_')

set.seed(10)
bm <- Seurat::FindMultiModalNeighbors(
  bm, reduction.list = list("pca", "apca"), 
  dims.list = list(1:30, 1:30), modality.weight.name = "RNA.weight"
)

set.seed(10)
bm <- Seurat::RunUMAP(bm, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

##########

plot1 <-Seurat::DimPlot(bm, reduction = "rna.umap",
                        group.by = "ct", label = TRUE,
                        repel = TRUE, label.size = 2.5,
                        cols = col_palette)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human BM (Abseq, RNA+ADT)\nRNA UMAP"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/abseq_bm97_rna-umap.png"),
                plot1, device = "png", width = 9, height = 5, units = "in")

plot1 <-Seurat::DimPlot(bm, reduction = "adt.umap",
                        group.by = "ct", label = TRUE,
                        repel = TRUE, label.size = 2.5,
                        cols = col_palette)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human BM (Abseq, RNA+ADT)\nADT UMAP"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/abseq_bm97_adt-umap.png"),
                plot1, device = "png", width = 9, height = 5, units = "in")

plot1 <-Seurat::DimPlot(bm, reduction = "wnn.umap",
                        group.by = "ct", label = TRUE,
                        repel = TRUE, label.size = 2.5,
                        cols = col_palette)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human BM (Abseq, RNA+ADT)\nWNN UMAP"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/abseq_bm97_wnn-umap.png"),
                plot1, device = "png", width = 9, height = 5, units = "in")

##########################

bm[["bothUMAP"]] <- NULL
bm[["bothTSNE"]] <- NULL
bm[["MOFA"]] <- NULL
bm[["MOFAUMAP"]] <- NULL
bm[["MOFATSNE"]] <- NULL
bm[["Projected"]] <- NULL
bm[["ProjectedMean"]] <- NULL
bm[["BOTH"]] <- NULL

save(bm, date_of_run, session_info,
     file = "../../../out/main/abseq_bm97_preprocessed.RData")

