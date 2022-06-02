rm(list=ls())

library(Seurat); library(Signac); library(slingshot)
source("slingshot_funcs.R")

load("../../../out/main/10x_greenleaf_tcca_RNA-ATAC.RData")
greenleaf_tmp <- greenleaf
load("../../../out/main/10x_greenleaf_preprocessed.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

greenleaf[["common_tccaUMAP"]] <- greenleaf_tmp[["common_tcca"]]

multiSVD_obj <- tiltedCCA:::tiltedCCA_decomposition(input_obj = multiSVD_obj,
                                                    verbose = 1,
                                                    bool_modality_1_full = F,
                                                    bool_modality_2_full = F)
multiSVD_obj <-  tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 1)
dimred_1 <-  tiltedCCA:::.get_tCCAobj(multiSVD_obj, apply_postDimred = F, what = "common_dimred")
multiSVD_obj <-  tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 2)
dimred_2 <-  tiltedCCA:::.get_tCCAobj(multiSVD_obj, apply_postDimred = F, what = "common_dimred")
dimred <- cbind(dimred_1, dimred_2)
colnames(dimred) <- paste0("common_tcca_", 1:ncol(dimred))
greenleaf[["common_tccaDimred"]] <- Seurat::CreateDimReducObject(dimred)

Seurat::DefaultAssay(greenleaf) <- "SCT"
set.seed(10)
greenleaf <- Seurat::FindNeighbors(greenleaf, dims = 1:50, reduction = "pca")
set.seed(10)
greenleaf <- Seurat::FindClusters(greenleaf, resolution = 1)
greenleaf$rna_clustering <- greenleaf$seurat_clusters

Seurat::DefaultAssay(greenleaf) <- "ATAC"
set.seed(10)
greenleaf <- Seurat::FindNeighbors(greenleaf, dims = 2:50, reduction = "lsi")
set.seed(10)
greenleaf <- Seurat::FindClusters(greenleaf, resolution = 1)
greenleaf$atac_clustering <- greenleaf$seurat_clusters

Seurat::DefaultAssay(greenleaf) <- "SCT"
set.seed(10)
greenleaf <- Seurat::FindNeighbors(greenleaf, dims = 1:49, 
                                   reduction = "consensusPCA",
                                   graph.name = c("consensus_nn", "consensus_snn"))
set.seed(10)
greenleaf <- Seurat::FindClusters(greenleaf, resolution = 1,
                                  graph.name = "consensus_snn")
greenleaf$consensus_clustering <- greenleaf$seurat_clusters

Seurat::DefaultAssay(greenleaf) <- "SCT"
p <- ncol(greenleaf[["common_tccaDimred"]]@cell.embeddings)
set.seed(10)
greenleaf <- Seurat::FindNeighbors(greenleaf, dims = 1:p, 
                                   reduction = "common_tccaDimred",
                                   graph.name = c("tcca_nn", "tcca_snn"))
set.seed(10)
greenleaf <- Seurat::FindClusters(greenleaf, resolution = 1,
                                  graph.name = "tcca_snn")
greenleaf$tcca_clustering <- greenleaf$seurat_clusters

table(greenleaf$rna_clustering, greenleaf$atac_clustering)
table(greenleaf$rna_clustering, greenleaf$consensus_clustering)
table(greenleaf$rna_clustering, greenleaf$tcca_clustering)
table(greenleaf$consensus_clustering, greenleaf$tcca_clustering)

##############################

plot1 <- Seurat::DimPlot(greenleaf, reduction = "umap",
                         group.by = "rna_clustering", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nRNA clustering"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_greenleaf_lineageAlignment_rna-clustering.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <- Seurat::DimPlot(greenleaf, reduction = "umap.atac",
                         group.by = "atac_clustering", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nATAC clustering"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_greenleaf_lineageAlignment_atac-clustering.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <- Seurat::DimPlot(greenleaf, reduction = "consensusUMAP",
                         group.by = "consensus_clustering", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nConsensus clustering"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_greenleaf_lineageAlignment_consensus-clustering.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <- Seurat::DimPlot(greenleaf, reduction = "common_tccaUMAP",
                         group.by = "tcca_clustering", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nTCCA clustering"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_greenleaf_lineageAlignment_tcca-clustering.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

##############################

dimred <- greenleaf[["pca"]]@cell.embeddings[,1:50]
lineage_order <- c("10", "7", "4", "0", "6", "5", "8", "9")
cluster_vec <- as.character(greenleaf$rna_clustering)
keep_idx <- which(cluster_vec %in% lineage_order)
dimred <- dimred[keep_idx,]
cluster_vec <- factor(cluster_vec[keep_idx], levels = sort(unique(cluster_vec[keep_idx])))
initial_fit <- .initial_curve_fit(cluster_vec = cluster_vec,
                                    dimred = dimred, 
                                    lineage_order = lineage_order)
pseudotime <- .extract_pseudotime(dimred = dimred,
                                  initial_fit = initial_fit,
                                  stretch = 2)
pseudotime <- rank(pseudotime)
pseudotime_full <- rep(NA, ncol(greenleaf))
names(pseudotime_full) <- colnames(greenleaf)
pseudotime_full[names(pseudotime)] <- pseudotime
greenleaf$rna_pseudotime <- pseudotime_full

plot1 <- Seurat::FeaturePlot(greenleaf, reduction = "umap",
                         features = "rna_pseudotime")
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nRNA pseudotime"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_greenleaf_lineageAlignment_rna-pseudotime.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

########################

dimred <- greenleaf[["lsi"]]@cell.embeddings[,2:50]
cluster_vec <- as.character(greenleaf$atac_clustering)
# cluster_vec[cluster_vec == "6"] <- "5"
lineage_order <- c("4", "0", "10", "6", "5", "8")
keep_idx <- which(cluster_vec %in% lineage_order)
dimred <- dimred[keep_idx,]
cluster_vec <- factor(cluster_vec[keep_idx], levels = sort(unique(cluster_vec[keep_idx])))
initial_fit <- .initial_curve_fit(cluster_vec = cluster_vec,
                                  dimred = dimred, 
                                  lineage_order = lineage_order)
pseudotime <- .extract_pseudotime(dimred = dimred,
                                  initial_fit = initial_fit,
                                  stretch = 2)
pseudotime <- rank(pseudotime)
pseudotime_full <- rep(NA, ncol(greenleaf))
names(pseudotime_full) <- colnames(greenleaf)
pseudotime_full[names(pseudotime)] <- pseudotime
greenleaf$atac_pseudotime <- pseudotime_full

plot1 <- Seurat::FeaturePlot(greenleaf, reduction = "umap.atac",
                             features = "atac_pseudotime")
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nATAC pseudotime"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_greenleaf_lineageAlignment_atac-pseudotime.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

################################################

# before doing consensus PCA and tiltedCCA, we first need to find all the cells that 
# will be participating in this analysis
rna_pseudotime <- greenleaf$rna_pseudotime
atac_pseudotime <- greenleaf$atac_pseudotime

cell_idx <- intersect(which(!is.na(rna_pseudotime)), which(!is.na(atac_pseudotime)))
cell_names_keep <- colnames(greenleaf)[cell_idx]

# replot consensus and tilted clustering
consensus_clustering <- greenleaf$consensus_clustering 
idx <- which(!names(consensus_clustering) %in% cell_names_keep)
consensus_clustering[idx] <- NA
greenleaf$consensus_clustering_subset <- consensus_clustering

plot1 <- Seurat::DimPlot(greenleaf, reduction = "consensusUMAP",
                         group.by = "consensus_clustering_subset", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nConsensus clustering"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_greenleaf_lineageAlignment_consensus-clustering_subset.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

tcca_clustering <- greenleaf$tcca_clustering 
idx <- which(!names(tcca_clustering) %in% cell_names_keep)
tcca_clustering[idx] <- NA
greenleaf$tcca_clustering_subset <- tcca_clustering

plot1 <- Seurat::DimPlot(greenleaf, reduction = "common_tccaUMAP",
                         group.by = "tcca_clustering_subset", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nTCCA clustering"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_greenleaf_lineageAlignment_tcca-clustering_subset.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

#########################

# now fit the new ordering for consensus PCA and tilted CCA
dimred <- greenleaf[["consensusPCA"]]@cell.embeddings[,1:49]
lineage_order <- c("9", "5", "0", "4", "2", "10")
cluster_vec <- as.character(greenleaf$consensus_clustering)
keep_idx <- which(cluster_vec %in% lineage_order)
dimred_subset1 <- dimred[keep_idx,]
cluster_vec <- factor(cluster_vec[keep_idx], levels = sort(unique(cluster_vec[keep_idx])))
initial_fit <- .initial_curve_fit(cluster_vec = cluster_vec,
                                  dimred = dimred_subset1, 
                                  lineage_order = lineage_order)
dimred_subset2 <- dimred[cell_names_keep,]
pseudotime <- .extract_pseudotime(dimred = dimred_subset2,
                                  initial_fit = initial_fit,
                                  stretch = 2)
pseudotime <- rank(pseudotime)
pseudotime_full <- rep(NA, ncol(greenleaf))
names(pseudotime_full) <- colnames(greenleaf)
pseudotime_full[names(pseudotime)] <- pseudotime
greenleaf$consensus_pseudotime <- pseudotime_full

plot1 <- Seurat::FeaturePlot(greenleaf, reduction = "consensusUMAP",
                             features = "consensus_pseudotime")
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nConsensus pseudotime"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_greenleaf_lineageAlignment_consensus-pseudotime.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

####

dimred <- greenleaf[["common_tccaDimred"]]@cell.embeddings
lineage_order <- c("12", "2", "8", "4", "1", "3", "10")
cluster_vec <- as.character(greenleaf$tcca_clustering)
keep_idx <- which(cluster_vec %in% lineage_order)
dimred_subset1 <- dimred[keep_idx,]
cluster_vec <- factor(cluster_vec[keep_idx], levels = sort(unique(cluster_vec[keep_idx])))
initial_fit <- .initial_curve_fit(cluster_vec = cluster_vec,
                                  dimred = dimred_subset1, 
                                  lineage_order = lineage_order)
dimred_subset2 <- dimred[cell_names_keep,]
pseudotime <- .extract_pseudotime(dimred = dimred_subset2,
                                  initial_fit = initial_fit,
                                  stretch = 2)
pseudotime <- rank(pseudotime)
pseudotime_full <- rep(NA, ncol(greenleaf))
names(pseudotime_full) <- colnames(greenleaf)
pseudotime_full[names(pseudotime)] <- pseudotime
greenleaf$tcca_pseudotime <- pseudotime_full

plot1 <- Seurat::FeaturePlot(greenleaf, reduction = "common_tccaUMAP",
                             features = "tcca_pseudotime")
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nTCCA pseudotime"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_greenleaf_lineageAlignment_tcca-pseudotime.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

save(greenleaf, date_of_run, session_info,
     file = "../../../out/main/10x_greenleaf_lineageAlignment_pseudotime.RData")

