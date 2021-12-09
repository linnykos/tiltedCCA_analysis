# see if the proteins I found are convincingly that distinct anyways
rm(list=ls())
load("../../../../out/Writeup14f/Writeup14f_citeseq_bm_dcca.RData")
dcca_decomp <- multiomicCCA::dcca_decomposition(dcca_res2)

n <- nrow(bm@meta.data)
svd_1 <- multiomicCCA:::.svd_truncated(dcca_decomp$common_mat_1, 
                                       K = rank_1, K_full_rank = F,
                                       rescale = F,
                                       mean_vec = F, sd_vec = F,
                                       symmetric = F)
svd_2 <- multiomicCCA:::.svd_truncated(dcca_decomp$common_mat_2, K = rank_2, K_full_rank = F,
                                       rescale = F,
                                       mean_vec = F, sd_vec = F,
                                       symmetric = F)
dimred_1 <- multiomicCCA:::.mult_mat_vec(svd_1$u, svd_1$d/svd_1$d[1]*sqrt(n))
l2_vec <- apply(dimred_1, 1, function(x){multiomicCCA:::.l2norm(x)})
dimred_1 <- multiomicCCA:::.mult_vec_mat(1/l2_vec, dimred_1)
dimred_2 <- multiomicCCA:::.mult_mat_vec(svd_2$u, svd_2$d/svd_2$d[1]*sqrt(n))
l2_vec <- apply(dimred_2, 1, function(x){multiomicCCA:::.l2norm(x)})
dimred_2 <- multiomicCCA:::.mult_vec_mat(1/l2_vec, dimred_2)

set.seed(10)
common_umap <- Seurat::RunUMAP(cbind(dimred_1, dimred_2), 
                               metric = "euclidean",
                               reduction.key = "umapCommon_")
rownames(common_umap@cell.embeddings) <- rownames(bm@meta.data)
bm[["dcca_common"]] <- Seurat::CreateDimReducObject(common_umap@cell.embeddings)

#######################################

svd_2 <- multiomicCCA:::.svd_truncated(dcca_decomp$distinct_mat_2, 
                                       K = rank_2, K_full_rank = F,
                                       rescale = F,
                                       mean_vec = F, sd_vec = F,
                                       symmetric = F)
dimred_2 <- multiomicCCA:::.mult_mat_vec(svd_2$u, svd_2$d)

set.seed(10)
distinct2_umap <- Seurat::RunUMAP(dimred_2, 
                                  metric = "cosine",
                                  reduction.key = "umapDistinct2_")
rownames(distinct2_umap@cell.embeddings) <- rownames(bm@meta.data)
bm[["dcca_distinct2"]] <- Seurat::CreateDimReducObject(distinct2_umap@cell.embeddings)

################################

anchor_name <- "rna.umap"
other_names <- c("adt.umap", "dcca_common", "dcca_distinct2")

for(umap_name in other_names){
  print(umap_name)
  u_mat1 <- bm[[anchor_name]]@cell.embeddings
  u_mat2 <- bm[[umap_name]]@cell.embeddings
  tmp <- svd(t(u_mat1) %*% u_mat2)
  rotation_mat <- tmp$u %*% t(tmp$v)
  tmp <- u_mat2 %*% t(rotation_mat)
  rownames(tmp) <- rownames(bm@meta.data)
  colnames(tmp) <- colnames(bm[[umap_name]]@cell.embeddings)
  bm[[umap_name]]@cell.embeddings <- tmp
}

#######################################

tmp <- dcca_decomp$common_mat_2
colnames(tmp) <- paste0("c-", colnames(tmp))
bm[["common_protein"]] <- Seurat::CreateAssayObject(counts = t(tmp))

tmp <- dcca_decomp$distinct_mat_2
colnames(tmp) <- paste0("d2-", colnames(tmp))
bm[["distinct_protein"]] <- Seurat::CreateAssayObject(counts = t(tmp))

reduction_vec <- c("adt.umap", "dcca_common", "dcca_distinct2")
assay_vec <- c("ADT", "common_protein", "distinct_protein")
prefix <- c("", "c-", "d2-")
protein_vec <- c("CD25", "CD197-CCR7", "CD57", "CD45RO")
main_vec <- c("(Protein)", 
              "(Tilted-CCA, Common)", 
              "(T-CCA, Distinct 2, Alignment)")
file_vec <- c("adt", "tiltedcca-common", "tiltedcca-distinct2")
for(i in 1:length(reduction_vec)){
  Seurat::DefaultAssay(bm) <- assay_vec[i]
  plot1 <- Seurat::FeaturePlot(bm, features = paste0(prefix[i], protein_vec),
                               reduction = reduction_vec[i])
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14f/Writeup14f_citeseq_bm_distinct_protein_", reduction_vec[i], ".png"),
                  plot1, device = "png", width = 9, height = 8, units = "in")
}

protein_vec <- c("CD14", "CD11c", "CD19", "CD3")
for(i in 1:length(reduction_vec)){
  Seurat::DefaultAssay(bm) <- assay_vec[i]
  plot1 <- Seurat::FeaturePlot(bm, features = paste0(prefix[i], protein_vec),
                               reduction = reduction_vec[i])
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14f/Writeup14f_citeseq_bm_common_protein_", reduction_vec[i], ".png"),
                  plot1, device = "png", width = 9, height = 8, units = "in")
}

protein_vec <- c("CD28", "CD3", "CD38", 
                 "CD4", "CD45RA", "CD8a",
                 "CD278-ICOS", "CD27", "CD11a")
for(i in 1:length(reduction_vec)){
  Seurat::DefaultAssay(bm) <- assay_vec[i]
  plot1 <- Seurat::FeaturePlot(bm, features = paste0(prefix[i], protein_vec),
                               reduction = reduction_vec[i])
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14f/Writeup14f_citeseq_bm_common_cd4cd8_", reduction_vec[i], ".png"),
                  plot1, device = "png", width = 12, height = 10, units = "in")
}

