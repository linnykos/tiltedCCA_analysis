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
bm <- Seurat::RunUMAP(bm, reduction = 'pca', dims = 1:50, assay = 'RNA',
                      reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
set.seed(10)
bm <- Seurat::RunUMAP(bm, 
                      reduction = 'apca', 
                      dims = 1:50, assay = 'ADT',
                      reduction.name = 'adt.umap', 
                      reduction.key = 'adtUMAP_')

set.seed(10)
bm <- Seurat::FindMultiModalNeighbors(
  bm, reduction.list = list("pca", "apca"), 
  dims.list = list(1:50, 1:50), modality.weight.name = "RNA.weight"
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
                                           dims_1 = 1:50, dims_2 = 1:50,
                                           dims_consensus = 1:50,
                                           verbose = 1)
consensus_dimred <- consensus_pca$dimred_consensus
colnames(consensus_dimred) <- paste0("consensusPCA_", 1:ncol(consensus_dimred))
bm[["consensusPCA"]] <- Seurat::CreateDimReducObject(consensus_dimred, 
                                                     assay = "RNA")

set.seed(10)
umap_res <- Seurat::RunUMAP(consensus_pca$dimred_consensus)
save(umap_res, "../../../out/main/tmp_consensuspca.RData")

umap_mat <- umap_res@cell.embeddings
rownames(umap_mat) <- colnames(bm)
bm[["consensusUMAP"]] <- Seurat::CreateDimReducObject(umap_mat, 
                                                      assay = "RNA",
                                                      reduction.name = 'consensus.umap', 
                                                      reduction.key = 'consensusUMAP_')

##########

plot1 <-Seurat::DimPlot(bm, reduction = "rna.umap",
                        group.by = "ct", label = TRUE,
                        repel = TRUE, label.size = 2.5,
                        cols = col_palette)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human BM (Abseq, RNA+ADT)\nRNA UMAP"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/abseq_bm97Ref_rna-umap.png"),
                plot1, device = "png", width = 11, height = 5, units = "in")

plot1 <-Seurat::DimPlot(bm, reduction = "adt.umap",
                        group.by = "ct", label = TRUE,
                        repel = TRUE, label.size = 2.5,
                        cols = col_palette)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human BM (Abseq, RNA+ADT)\nADT UMAP"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/abseq_bm97Ref_adt-umap.png"),
                plot1, device = "png", width = 11, height = 5, units = "in")

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
     file = "../../../out/main/abseq_bm97Ref_preprocessed_tmp.RData")

