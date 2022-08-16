rm(list=ls())
library(Seurat)
library(tiltedCCA)

load("../../../out/main/citeseq_bm25_preprocessed.RData")
source("bm_25antibody_colorPalette.R")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

# Seurat::DefaultAssay(bm) <- "RNA"
# bm <- Seurat::ScaleData(bm) 
# bm <- Seurat::RunPCA(bm, verbose = F)
# Seurat::DefaultAssay(bm) <- "ADT"
# bm <- Seurat::ScaleData(bm)
# bm <- Seurat::RunPCA(bm, reduction.name = 'apca', verbose = F)
# 
# pca_rna <- bm[["pca"]]@cell.embeddings[,1:30]
# pca_adt <- bm[["apca"]]@cell.embeddings[,1:18]
# n <- nrow(pca_rna)
# sing_rna <- apply(pca_rna, 2, function(x){sqrt(sum(x^2))})
# sing_adt <- apply(pca_adt, 2, function(x){sqrt(sum(x^2))})
# pca_rna <- tiltedCCA:::.mult_mat_vec(pca_rna, sqrt(n)/sing_rna)
# pca_adt <- tiltedCCA:::.mult_mat_vec(pca_adt, sqrt(n)/sing_adt)
# 
# consensus_mat <- cbind(pca_rna, pca_adt)
# set.seed(10)
# umap_res <- Seurat::RunUMAP(consensus_mat)
# 
# umap_mat <- umap_res@cell.embeddings
# rownames(umap_mat) <- colnames(bm)
# bm[["consensusPCA"]] <- Seurat::CreateDimReducObject(umap_mat)

##############################

# consensus pca
Seurat::DefaultAssay(bm) <- "RNA"
mat_1 <- Matrix::t(bm[["RNA"]]@data[Seurat::VariableFeatures(object = bm),])
Seurat::DefaultAssay(bm) <- "ADT"
mat_2 <- Matrix::t(bm[["ADT"]]@data[Seurat::VariableFeatures(object = bm),])

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
consensus_pca <- tiltedCCA:::consensus_pca(mat_1 = mat_1b, mat_2 = mat_2b,
                                           dims_1 = 1:30, dims_2 = 1:18,
                                           dims_consensus = 1:30,
                                           center_1 = T, center_2 = T,
                                           recenter_1 = F, recenter_2 = F,
                                           rescale_1 = F, rescale_2 = F,
                                           scale_1 = T, scale_2 = T,
                                           verbose = 1)
consensus_dimred <- consensus_pca$dimred_consensus
colnames(consensus_dimred) <- paste0("consensusPCA_", 1:ncol(consensus_dimred))
bm[["consensusPCA"]] <- Seurat::CreateDimReducObject(consensus_dimred, 
                                                     assay = "RNA")

set.seed(10)
umap_res <- Seurat::RunUMAP(consensus_pca$dimred_consensus)
umap_mat <- umap_res@cell.embeddings
rownames(umap_mat) <- colnames(bm)
colnames(umap_mat) <- paste0("consensusUMAP_", 1:ncol(umap_mat))
bm[["consensusUMAP"]] <- Seurat::CreateDimReducObject(umap_mat, 
                                                      assay = "RNA")

plot1 <- Seurat::DimPlot(bm, reduction = "consensusUMAP",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         cols = col_palette)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human BM (CITE-Seq, RNA+ADT)\nConsensus PCA"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_bm25_consensusPCA.png"),
                plot1, device = "png", width = 5.5, height = 5, units = "in")

plot1 <- Seurat::DimPlot(bm, reduction = "consensusUMAP",
                         group.by = "celltype.l2", 
                         cols = col_palette)
plot1 <- plot1 + Seurat::NoLegend() + Seurat::NoAxes()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_bm25_consensusPCA_cleaned.png"),
                plot1, device = "png", width = 3, height = 3, units = "in",
                dpi = 500)

############################################################

dimred <- consensus_pca$dimred_1
for(j in 1:ncol(dimred)){
  tmp_df <- data.frame(cbind(dimred[,j], consensus_pca$dimred_consensus))
  colnames(tmp_df)[1] <- "y"
  lm_res <- stats::lm(y ~ . - 1, data = tmp_df)
  dimred[,j] <- stats::residuals(lm_res)
}

set.seed(10)
umap_res <- Seurat::RunUMAP(dimred)
umap_mat <- umap_res@cell.embeddings
rownames(umap_mat) <- colnames(bm)
colnames(umap_mat) <- paste0("distinctRNA_", 1:ncol(umap_mat))
bm[["distinctRNA"]] <- Seurat::CreateDimReducObject(umap_mat, 
                                                    assay = "RNA")

plot1 <- Seurat::DimPlot(bm, reduction = "distinctRNA",
                         group.by = "celltype.l2", 
                         cols = col_palette)
plot1 <- plot1 + Seurat::NoLegend() + Seurat::NoAxes()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_bm25_consensusPCA-distinctRNA_cleaned.png"),
                plot1, device = "png", width = 3, height = 3, units = "in",
                dpi = 500)

#####################3

dimred <- consensus_pca$dimred_2
for(j in 1:ncol(dimred)){
  tmp_df <- data.frame(cbind(dimred[,j], consensus_pca$dimred_consensus))
  colnames(tmp_df)[1] <- "y"
  lm_res <- stats::lm(y ~ . - 1, data = tmp_df)
  dimred[,j] <- stats::residuals(lm_res)
}

set.seed(10)
umap_res <- Seurat::RunUMAP(dimred)
umap_mat <- umap_res@cell.embeddings
rownames(umap_mat) <- colnames(bm)
colnames(umap_mat) <- paste0("distinctADT_", 1:ncol(umap_mat))
bm[["distinctADT"]] <- Seurat::CreateDimReducObject(umap_mat, 
                                                    assay = "RNA")

plot1 <- Seurat::DimPlot(bm, reduction = "distinctADT",
                         group.by = "celltype.l2", 
                         cols = col_palette)
plot1 <- plot1 + Seurat::NoLegend() + Seurat::NoAxes()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_bm25_consensusPCA-distinctADT_cleaned.png"),
                plot1, device = "png", width = 3, height = 3, units = "in",
                dpi = 500)
