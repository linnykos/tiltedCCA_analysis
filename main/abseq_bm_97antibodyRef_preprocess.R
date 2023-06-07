rm(list=ls())
source("bm_97antibodyRef_colorPalette.R")

library(Seurat)
library(dplyr)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

bm <- readRDS("~/nzhanglab/data/Triana_nature_bonemarrow/Healthy.rds")
bm$ct <- factor(bm$ct, levels = sort(levels(bm$ct)))

##################

DefaultAssay(bm) <- "RNA"
bm[["RNA"]]@var.features <- rownames(bm)
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

###############

# consensus pca
Seurat::DefaultAssay(bm) <- "RNA"
mat_1 <- Matrix::t(bm[["RNA"]]@scale.data[Seurat::VariableFeatures(object = bm),])
Seurat::DefaultAssay(bm) <- "AB"
mat_2 <- Matrix::t(bm[["AB"]]@scale.data)

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

consensus_pca <- tiltedCCA:::consensus_pca(mat_1 = mat_1b, mat_2 = mat_2b,
                                           dims_1 = 1:30, dims_2 = 1:30,
                                           dims_consensus = 1:20,
                                           verbose = 1)
consensus_dimred <- consensus_pca$dimred_consensus
colnames(consensus_dimred) <- paste0("consensusPCA_", 1:ncol(consensus_dimred))
bm[["consensusPCA"]] <- Seurat::CreateDimReducObject(consensus_dimred, 
                                                     assay = "RNA")

set.seed(10)
umap_res <- Seurat::RunUMAP(consensus_pca$dimred_consensus)
# save(umap_res, file = "../../../out/main/tmp_consensuspca.RData")

umap_mat <- umap_res@cell.embeddings
rownames(umap_mat) <- colnames(bm)
colnames(umap_mat) <- paste0("consensusUMAP_", 1:ncol(umap_mat))
bm[["consensusUMAP"]] <- Seurat::CreateDimReducObject(umap_mat, 
                                                      assay = "RNA")

##########

plot1 <-Seurat::DimPlot(bm, reduction = "rna.umap",
                        group.by = "ct", label = TRUE,
                        repel = TRUE, label.size = 2.5,
                        cols = col_palette)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human BM (Abseq, RNA+ADT)\nRNA UMAP"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/abseq_bm97Ref_rna-umap.png"),
                plot1, device = "png", width = 11, height = 5, units = "in")

plot1 <- Seurat::DimPlot(bm, reduction = "rna.umap",
                         group.by = "ct", 
                         cols = col_palette)
plot1 <- plot1 + Seurat::NoLegend() + Seurat::NoAxes()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/abseq_bm97Ref_rna-umap_cleaned.png"),
                plot1, device = "png", width = 3, height = 3, units = "in",
                dpi = 500)

plot1 <-Seurat::DimPlot(bm, reduction = "adt.umap",
                        group.by = "ct", label = TRUE,
                        repel = TRUE, label.size = 2.5,
                        cols = col_palette)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human BM (Abseq, RNA+ADT)\nADT UMAP"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/abseq_bm97Ref_adt-umap.png"),
                plot1, device = "png", width = 11, height = 5, units = "in")

plot1 <- Seurat::DimPlot(bm, reduction = "adt.umap",
                         group.by = "ct", 
                         cols = col_palette)
plot1 <- plot1 + Seurat::NoLegend() + Seurat::NoAxes()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/abseq_bm97Ref_adt-umap_cleaned.png"),
                plot1, device = "png", width = 3, height = 3, units = "in",
                dpi = 500)

plot1 <-Seurat::DimPlot(bm, reduction = "wnn.umap",
                        group.by = "ct", label = TRUE,
                        repel = TRUE, label.size = 2.5,
                        cols = col_palette)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human BM (Abseq, RNA+ADT)\nWNN UMAP"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/abseq_bm97Ref_wnn-umap.png"),
                plot1, device = "png", width = 11, height = 5, units = "in")

plot1 <-Seurat::DimPlot(bm, reduction = "consensusUMAP",
                        group.by = "ct", label = TRUE,
                        repel = TRUE, label.size = 2.5,
                        cols = col_palette)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human BM (Abseq, RNA+ADT)\nConsensus PCA's UMAP"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/abseq_bm97Ref_consensusPCA-umap.png"),
                plot1, device = "png", width = 11, height = 5, units = "in")

plot1 <- Seurat::DimPlot(bm, reduction = "consensusUMAP",
                         group.by = "ct", 
                         cols = col_palette)
plot1 <- plot1 + Seurat::NoLegend() + Seurat::NoAxes()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/abseq_bm97Ref_consensusPCA-umap_cleaned.png"),
                plot1, device = "png", width = 3, height = 3, units = "in",
                dpi = 500)

##########################

bm[["UMAPni"]] <- NULL
bm[["TSNEni"]] <- NULL
bm[["tsne"]] <- NULL
bm[["MOFA"]] <- NULL
bm[["MOFAUMAP"]] <- NULL
bm[["MOFATSNE"]] <- NULL
bm[["Projected"]] <- NULL
bm[["ProjectedMean"]] <- NULL
bm[["BOTH"]] <- NULL
bm[["integrated"]] <- NULL

Seurat::DefaultAssay(bm) <- "RNA"
save(bm, date_of_run, session_info,
     file = "../../../out/main/abseq_bm97Ref_preprocessed.RData")

##########################

source("bm_97antibodyRef_colorPalette.R")

bm[["adt.umap"]]@cell.embeddings <- tiltedCCA:::.rotate_matrix(
  source_mat = bm[["rna.umap"]]@cell.embeddings,
  target_mat = bm[["adt.umap"]]@cell.embeddings
)
plot1 <- Seurat::DimPlot(bm, reduction = "adt.umap",
                         group.by = "ct", 
                         cols = col_palette, pt.size = 1.5)
plot1 <- plot1 + Seurat::NoLegend() 
plot1 <- plot1 + ggplot2::ggtitle("") + ggplot2::ylab("AB UMAP 2") + ggplot2::xlab("AB UMAP 1")
ggplot2::ggsave(filename = paste0("../../../out/figures/main/abseq_bm97Ref_adt-umap_slides.png"),
                plot1, device = "png", width = 4, height = 4, units = "in")

plot1 <- Seurat::DimPlot(bm, reduction = "rna.umap",
                         group.by = "ct", 
                         cols = col_palette, pt.size = 1.5)
plot1 <- plot1 + Seurat::NoLegend() 
plot1 <- plot1 + ggplot2::ggtitle("") + ggplot2::ylab("RNA UMAP 2") + ggplot2::xlab("RNA UMAP 1")
ggplot2::ggsave(filename = paste0("../../../out/figures/main/abseq_bm97Ref_rna-umap_slides.png"),
                plot1, device = "png", width = 4, height = 4, units = "in")

gray_vec <- rep("gray50", length(unique(bm$ct)))
names(gray_vec) <- unique(bm$ct)
plot1 <- Seurat::DimPlot(bm, reduction = "adt.umap", group.by = "ct", cols = gray_vec, pt.size = 0.1)
plot1 <- plot1 + Seurat::NoLegend() 
plot1 <- plot1 + ggplot2::ggtitle("") + ggplot2::ylab("AB UMAP 2") + ggplot2::xlab("AB UMAP 1")
ggplot2::ggsave(filename = paste0("../../../out/figures/main/abseq_bm97Ref_adt-umap_slides2.png"),
                plot1, device = "png", width = 4, height = 4, units = "in")

plot1 <- Seurat::DimPlot(bm, reduction = "rna.umap", group.by = "ct", cols = gray_vec, pt.size = 0.1)
plot1 <- plot1 + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("") + ggplot2::ylab("RNA UMAP 2") + ggplot2::xlab("RNA UMAP 1")
ggplot2::ggsave(filename = paste0("../../../out/figures/main/abseq_bm97Ref_rna-umap_slides2.png"),
                plot1, device = "png", width = 4, height = 4, units = "in")


bm[["consensusUMAP"]]@cell.embeddings <- tiltedCCA:::.rotate_matrix(
  source_mat = bm[["rna.umap"]]@cell.embeddings,
  target_mat = bm[["consensusUMAP"]]@cell.embeddings
)
plot1 <- Seurat::DimPlot(bm, reduction = "consensusUMAP",
                         group.by = "ct", 
                         cols = col_palette, pt.size = 1.5)
plot1 <- plot1 + Seurat::NoLegend() 
plot1 <- plot1 + ggplot2::ggtitle("") + ggplot2::ylab("ConsensusPCA UMAP 2") + ggplot2::xlab("ConsensusPCA UMAP 1")
ggplot2::ggsave(filename = paste0("../../../out/figures/main/abseq_bm97Ref_consensusPCA-umap_slides.png"),
                plot1, device = "png", width = 4, height = 4, units = "in")


bm[["consensusUMAP"]]@cell.embeddings <- tiltedCCA:::.rotate_matrix(
  source_mat = bm[["rna.umap"]]@cell.embeddings,
  target_mat = bm[["consensusUMAP"]]@cell.embeddings
)
plot1 <- Seurat::DimPlot(bm, reduction = "consensusUMAP",
                         group.by = "ct", 
                         cols = col_palette, pt.size = 0.05)
plot1 <- plot1 + Seurat::NoLegend() 
plot1 <- plot1 + ggplot2::ggtitle("") + ggplot2::ylab("ConsensusPCA UMAP 2") + ggplot2::xlab("ConsensusPCA UMAP 1")
ggplot2::ggsave(filename = paste0("../../../out/figures/main/abseq_bm97Ref_consensusPCA-umap_slides-small.png"),
                plot1, device = "png", width = 4, height = 4, units = "in")


