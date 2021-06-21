rm(list=ls())
load("../../../../out/Writeup11/Writeup11_citeseq_bm25_preprocessed.RData")
# load("../../../../out/Writeup13/Writeup13_citeseq_bm25_dcca.RData")

library(Seurat)
library(multiomicCCA)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

mat_1 <- t(bm[["RNA"]]@scale.data)
mat_2 <- t(bm[["ADT"]]@scale.data)

range(table(bm@meta.data$celltype.l2))

set.seed(10)
rank_1 <- 30; rank_2 <- 18
dcca_res <- multiomicCCA::dcca_factor(mat_1, mat_2, rank_1 = rank_1, rank_2 = rank_2, 
                                      meta_clustering = NA, num_neigh = 100,
                                      apply_shrinkage = F, fix_distinct_perc = F, 
                                      verbose = T) 
dcca_decomp <- multiomicCCA::dcca_decomposition(dcca_res, rank_c = min(rank_1, rank_2), 
                                                verbose = T)

## see the supplement of "Comprehensive Integration of Single-Cell Data", https://www.cell.com/cms/10.1016/j.cell.2019.05.031/attachment/1620c67d-66f0-476d-92e2-846b9f893846/mmc2.pdf for the marker genes

marker_genes <- c("AVP", "LMO4", "PF4", "BLVRB", "MME", "DERL3", "CLEC9A", "CDC1", "MPO", "AZU1", "CD14", "FCGR3A", "VREB3", "MS4A1", "CD79A", "IGKC", "PF4", "XCL1", "CD8A", "CD4", "SH2D1A")
marker_idx <- which(marker_genes %in% colnames(mat_1))
marker_genes[marker_idx]
marker_genes[-marker_idx]
marker_genes1 <- marker_genes[marker_idx]

membership_vec <- as.factor(bm@meta.data$celltype.l2)
set.seed(10)
zz <- multiomicCCA::plot_embeddings(dcca_decomp, membership_vec, data_1 = T, data_2 = F, 
                                    add_noise = T, pca = F, only_embedding = T)

bm_rna_com <- Seurat::CreateSeuratObject(counts = t(dcca_decomp$common_mat_1))
bm_rna_com[["rna_common"]] <- Seurat::CreateDimReducObject(embedding = zz[[1]], key = "commonUMAP_")

bm_rna_dis <- Seurat::CreateSeuratObject(counts = t(dcca_decomp$distinct_mat_1))
bm_rna_dis[["rna_distinct"]] <- Seurat::CreateDimReducObject(embedding = zz[[2]], key = "distinctUMAP_")

bm[["rna_everything"]] <- Seurat::CreateDimReducObject(embedding = zz[[3]], key = "everythingUMAP_")

## see https://stackoverflow.com/questions/48875135/save-multiple-plots-from-ggplot2-using-a-for-loop-by-side-by-side
for(i in 1:length(marker_genes1)){
  print(i)
  p1 <- Seurat::FeaturePlot(bm_rna_com, features = paste0("rna_", marker_genes1[i]), 
                            reduction = 'rna_common', ncol = 1)
  p2 <- Seurat::FeaturePlot(bm_rna_dis, features = paste0("rna_", marker_genes1[i]), 
                            reduction = 'rna_distinct', ncol = 1)
  p3 <- Seurat::FeaturePlot(bm, features = paste0("rna_", marker_genes1[i]), 
                            reduction = 'rna_everything', ncol = 1)
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14/Writeup14_citeseq_bm25_dcca_rna_marker", marker_genes1[i], ".png"),
                  gridExtra::arrangeGrob(p1,p2,p3,nrow = 1, ncol = 3), device = "png", width = 9, height = 3, units = "in")
}

main_vec <- c("common", "distinct")
reduction_vec <- c("rna_common", "rna_distinct")
for(i in 1:2){
  plot1 <- Seurat::DimPlot(bm, reduction = reduction_vec[i], group.by = 'celltype.l2', label = TRUE, 
                           repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
  plot1 <- plot1 + ggplot2::ggtitle(paste0("RNA ", main_vec[i], "  view (25, D-CCA)"))
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14/Writeup14_citeseq_bm25_dcca_rna_", main_vec[i], "_umap.png"),
                  plot1, device = "png", width = 5, height = 5, units = "in")
}

#################
# 
# membership_vec <- as.factor(bm@meta.data$celltype.l2)
# set.seed(10)
# zz <- multiomicCCA::plot_embeddings(dcca_decomp, membership_vec, data_1 = F, data_2 = T, 
#                                     add_noise = T, pca = F, only_embedding = T)
# bm[["adt_common"]] <- Seurat::CreateDimReducObject(embedding = zz[[1]], key = "adtcommonUMAP_")
# bm[["adt_distinct"]] <- Seurat::CreateDimReducObject(embedding = zz[[2]], key = "adtdistinctUMAP_")
# bm[["adt_everything"]] <- Seurat::CreateDimReducObject(embedding = zz[[3]], key = "adteverythingUMAP_")
# 
# main_vec <- c("common", "distinct", "everything")
# reduction_vec <- c("adt_common", "adt_distinct", "adt_everything")
# for(i in 1:3){
#   plot1 <- Seurat::DimPlot(bm, reduction = reduction_vec[i], group.by = 'celltype.l2', label = TRUE, 
#                            repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
#   plot1 <- plot1 + ggplot2::ggtitle(paste0("Protein ", main_vec[i], "  view (25, D-CCA)"))
#   ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14/Writeup14_citeseq_bm25_dcca_protein_", main_vec[i], "_umap.png"),
#                   plot1, device = "png", width = 5, height = 5, units = "in")
# }
# 
# for(i in 1:ncol(mat_2)){
#   print(i)
#   p1 <- Seurat::FeaturePlot(bm, features = paste0("adt_", colnames(mat_2)[i]), 
#                             reduction = 'adt_common', ncol = 1)
#   p2 <- Seurat::FeaturePlot(bm, features = paste0("adt_", colnames(mat_2)[i]), 
#                             reduction = 'adt_distinct', ncol = 1)
#   ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14/Writeup14_citeseq_bm25_dcca_adt_marker", colnames(mat_2)[i], ".png"),
#                   gridExtra::arrangeGrob(p1,p2, nrow = 1, ncol = 2), device = "png", width = 10, height = 5, units = "in")
# }
# 
