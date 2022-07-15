rm(list=ls())
library(Seurat)
library(Signac)
library(ArchR)
source("reik_colorPalette.R")

# reik <- readRDS("~/nzhanglab/data/GSE205117_Reik_mouseembryo/ftp1.babraham.ac.uk/data/processed/rna/seurat.rds")
# metadata <- read.csv("~/nzhanglab/data/GSE205117_Reik_mouseembryo/ftp1.babraham.ac.uk/sample_metadata.txt",
#                      sep = "\t")
# 
# set.seed(10)
# date_of_run <- Sys.time()
# session_info <- devtools::session_info()
# 
# n <- ncol(reik)
# mapping_vec <- rep(NA, n)
# for(i in 1:n){
#   if(i %% floor(n/10) == 0) cat('*')
#   tmp <- intersect(which(metadata$barcode == reik$barcode[i]),
#                    which(metadata$sample == reik$sample[i]))
#   
#   if(length(tmp) > 1) stop()
#   if(length(tmp) > 0) mapping_vec[i] <- tmp
# }
# 
# reik$stage <- metadata$stage[mapping_vec]
# reik$genotype <- metadata$genotype[mapping_vec]
# reik$pass_rnaQC <- metadata$pass_rnaQC[mapping_vec]
# reik$doublet_score <- metadata$doublet_score[mapping_vec]
# reik$doublet_call <- metadata$doublet_call[mapping_vec]
# reik$celltype <- metadata$celltype[mapping_vec]
# reik$celltype.score <- metadata$celltype.score[mapping_vec]
# reik$TSSEnrichment_atac <- metadata$TSSEnrichment_atac[mapping_vec]
# reik$ReadsInTSS_atac <- metadata$ReadsInTSS_atac[mapping_vec]
# reik$PromoterRatio_atac <- metadata$PromoterRatio_atac[mapping_vec]
# reik$NucleosomeRatio_atac <- metadata$NucleosomeRatio_atac[mapping_vec]
# reik$BlacklistRatio_atac <- metadata$BlacklistRatio_atac[mapping_vec]
# reik$pass_atacQC <- metadata$pass_atacQC[mapping_vec]
# reik$celltype.predicted <- metadata$celltype.predicted[mapping_vec]
# 
# set.seed(10)
# Seurat::DefaultAssay(reik) <- "RNA"
# reik <- Seurat::NormalizeData(reik)
# reik <- Seurat::FindVariableFeatures(reik)
# reik <- Seurat::ScaleData(reik)
# reik <- Seurat::RunPCA(reik, verbose = F)
# set.seed(10)
# reik <- Seurat::RunUMAP(reik, dims = 1:50)
# 
# save(reik, date_of_run, session_info,
#      file = "../../../out/main/10x_reik_preprocessed.RData")
# 
# ############
# 
# reik_atac <- ArchR::loadArchRProject(path = "~/nzhanglab/data/GSE205117_Reik_mouseembryo/ftp1.babraham.ac.uk/data/processed/atac/archR")
# reik_atac <- ArchR::addPeakMatrix(reik_atac, binarize = FALSE, force = TRUE) 
# ArchR::getAvailableMatrices(reik_atac) 
# 
# # This is a SummarizedExperiment object
# se_obj <- ArchR::getMatrixFromProject(ArchRProj = reik_atac,
#                                       useMatrix = "PeakMatrix") # this takes about 25 minutes
# 
# peak_granges <- se_obj@colData
# n <- ncol(reik)
# mapping_vec <- rep(NA, n)
# for(i in 1:n){
#   if(i %% floor(n/10) == 0) cat('*')
#   tmp <- intersect(which(peak_granges$barcode == reik$barcode[i]),
#                    which(peak_granges$sample == reik$sample[i]))
#   
#   if(length(tmp) > 1) stop()
#   if(length(tmp) > 0) mapping_vec[i] <- tmp 
# }
# 
# keep_vec <- rep(1, ncol(reik))
# keep_vec[which(is.na(mapping_vec))] <- 0
# reik$keep <- keep_vec
# reik <- subset(reik, keep == 1)
# 
# mapping_vec2 <- mapping_vec[which(!is.na(mapping_vec))]
# mat <- se_obj@assays@data$PeakMatrix[,mapping_vec2]
# chromosome_vec <- as.character(se_obj@rowRanges@seqnames)
# range_vec <- as.character(se_obj@rowRanges@ranges)
# p <- length(chromosome_vec)
# row_vec <- sapply(1:p, function(j){
#   paste0(chromosome_vec[j], "_", range_vec[j])
# })
# rownames(mat) <- row_vec
# reik[["ATAC"]] <- Seurat::CreateAssayObject(counts = mat)
# 
# Seurat::DefaultAssay(reik) <- "ATAC"
# set.seed(10)
# reik <- Signac::RunTFIDF(reik)
# reik <-  Signac::FindTopFeatures(reik, min.cutoff="q10")
# reik <-  Signac::RunSVD(reik)  
# 
# set.seed(10)
# reik <- Seurat::RunUMAP(reik, reduction="lsi", 
#                         dims=2:50, reduction.name="umap.atac", 
#                         reduction.key="atacUMAP_")
# 
# save(reik, peak_granges, date_of_run, session_info,
#      file = "../../../out/main/10x_reik_preprocessed.RData")
# 
# ###############
# 
# set.seed(10)
# reik <- Seurat::FindMultiModalNeighbors(reik, reduction.list = list("pca", "lsi"), 
#                                         dims.list = list(1:50, 2:50))
# set.seed(10)
# reik <- Seurat::RunUMAP(reik, nn.name = "weighted.nn", reduction.name = "umap.wnn", 
#                         reduction.key = "wnnUMAP_")
# 
# ###########

Seurat::DefaultAssay(reik) <- "RNA"
mat_1 <- Matrix::t(reik[["RNA"]]@data[Seurat::VariableFeatures(object = reik),])
Seurat::DefaultAssay(reik) <- "ATAC"
mat_2 <- Matrix::t(reik[["ATAC"]]@data[Seurat::VariableFeatures(object = reik),])

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
                                           dims_1 = 1:50, dims_2 = 2:50,
                                           dims_consensus = 1:49,
                                           center_1 = T, center_2 = F,
                                           recenter_1 = F, recenter_2 = T,
                                           rescale_1 = F, rescale_2 = T,
                                           scale_1 = T, scale_2 = F,
                                           verbose = 1)
consensus_dimred <- consensus_pca$dimred_consensus
colnames(consensus_dimred) <- paste0("consensusPCA_", 1:ncol(consensus_dimred))
reik[["consensusPCA"]] <- Seurat::CreateDimReducObject(consensus_dimred, 
                                                       assay = "RNA")

set.seed(10)
umap_res <- Seurat::RunUMAP(consensus_pca$dimred_consensus)
umap_mat <- umap_res@cell.embeddings
rownames(umap_mat) <- colnames(reik)
colnames(umap_mat) <- paste0("consensusUMAP_", 1:ncol(umap_mat))
reik[["consensusUMAP"]] <- Seurat::CreateDimReducObject(umap_mat, 
                                                        assay = "RNA")

save(reik, peak_granges, date_of_run, session_info,
     file = "../../../out/main/10x_reik_preprocessed.RData")

##############################

plot1 <- Seurat::DimPlot(reik, reduction = "umap",
                         group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         cols = col_palette)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Mouse embryo (10x, RNA+ATAC)\nRNA UMAP"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_reik_rna-umap_celltype.png"),
                plot1, device = "png", width = 10, height = 5, units = "in")

plot1 <- Seurat::DimPlot(reik, reduction = "umap",
                         group.by = "orig.ident", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Mouse embryo (10x, RNA+ATAC)\nRNA UMAP"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_reik_rna-umap_orig-ident.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

##

plot1 <- Seurat::DimPlot(reik, reduction = "umap.atac",
                         group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         cols = col_palette)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Mouse embryo (10x, RNA+ATAC)\nATAC UMAP"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_reik_atac-umap_celltype.png"),
                plot1, device = "png", width = 10, height = 5, units = "in")

plot1 <- Seurat::DimPlot(reik, reduction = "umap.atac",
                         group.by = "orig.ident", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Mouse embryo (10x, RNA+ATAC)\nATAC UMAP"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_reik_atac-umap_orig-ident.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

##

plot1 <- Seurat::DimPlot(reik, reduction = "umap.wnn",
                         group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         cols = col_palette)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Mouse embryo (10x, RNA+ATAC)\nWNN UMAP"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_reik_wnn-umap_celltype.png"),
                plot1, device = "png", width = 10, height = 5, units = "in")

plot1 <- Seurat::DimPlot(reik, reduction = "umap.wnn",
                         group.by = "orig.ident", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Mouse embryo (10x, RNA+ATAC)\nWNN UMAP"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_reik_wnn-umap_orig-ident.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

##

plot1 <- Seurat::DimPlot(reik, reduction = "consensusUMAP",
                         group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         cols = col_palette)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Mouse embryo (10x, RNA+ATAC)\nWNN UMAP"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_reik_consensuspca-umap_celltype.png"),
                plot1, device = "png", width = 10, height = 5, units = "in")

plot1 <- Seurat::DimPlot(reik, reduction = "consensusUMAP",
                         group.by = "orig.ident", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Mouse embryo (10x, RNA+ATAC)\nWNN UMAP"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_reik_consensuspca-umap_orig-ident.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")










