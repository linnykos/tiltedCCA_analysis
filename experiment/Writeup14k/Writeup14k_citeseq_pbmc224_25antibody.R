rm(list=ls())
load("../../../../out/Writeup14k/Writeup14k_citeseq_bm25_preprocessed.RData")
library(Seurat)
load("../../../../out/Writeup14j/Writeup14j_citeseq_pbmc224_dcca2.RData")
load("../../../../out/Writeup14i/Writeup14i_citeseq_pbmc224_dcca_variable_full_de.RData")
dcca_res <- dcca_res2
source("../Writeup14f/gene_exploration.R")

dcca_decomp <- multiomicCCA::dcca_decomposition(dcca_res)

summary_mat <- compute_variable_summary(mat = mat_2b, 
                                        common_mat = dcca_decomp$common_mat_2,
                                        metacell_clustering = NA)
colnames(summary_mat) <- c("r_squared", "p_value")
for(i in 1:nrow(summary_mat)){
  gene_name <- rownames(summary_mat)[i]
  
  max_celltype_idx <- max(combn_mat)
  celltype_pval <- sapply(1:max_celltype_idx, function(celltype){
    idx <- which(combn_mat == celltype, arr.ind = T)[,2]
    stats::quantile(sapply(idx, function(j){
      idx <- which(rownames(de_list[[j]]) == gene_name)
      if(length(idx) == 0) return(1)
      de_list[[j]][idx, "p_val"]
    }), probs = 0.75)
  })
  
  summary_mat[i,"p_value"] <- min(celltype_pval)
}
summary_mat[,"p_value"] <- -log10(summary_mat[,"p_value"])

################

set.seed(10)
pbmc2 <- pbmc
pbmc2[["adt.umap"]] <- NULL
pbmc2[["apca"]] <- NULL

# compute new ADT PCA
antibody_vec <- rownames(bm[["ADT"]])
print("Computing ADT PCA")
Seurat::DefaultAssay(pbmc2) <- "ADT"
antibody_vec[which(!antibody_vec %in% rownames(pbmc2[["ADT"]]@counts))]
sort(rownames(pbmc2[["ADT"]]@counts))
antibody_vec2 <- antibody_vec[which(antibody_vec %in% rownames(pbmc2[["ADT"]]@counts))]
antibody_vec2 <- c(antibody_vec2, "CD11a/CD18", "CD127", "CD278", "CD3-1", "CD3-2",
                   "CD38-1", "CD38-2", "CD4-1", "CD4-2", "CD56-1", "CD56-2", "HLA-DR")

df <- as.data.frame(summary_mat)
df$Labeling <- rep(0, nrow(df))
df$Labeling[which(rownames(df) %in% antibody_vec2)] <- 1
df$Labeling <- as.factor(df$Labeling)
df$Name <- rownames(df)
p1 <- ggplot2::ggplot(df, ggplot2::aes(x = p_value, y = r_squared))
p1 <- p1 + ggplot2::geom_point(ggplot2::aes(color = Labeling))
p1 <- p1 + ggplot2::scale_color_manual(values = c("black", "red"))
p1 <- p1 + ggplot2::scale_x_log10()
p1 <- p1 + ggplot2::ylim(0, 1)
p1 <- p1 + ggrepel::geom_text_repel(data = subset(df, Labeling == 1), 
                                    ggplot2::aes(label = Name, color = Labeling),
                                    max.overlaps = Inf, 
                                    size = 2)
p1 <- p1 + Seurat::NoLegend()
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14k/Writeup14k_citeseq_pbmc224_adt_exploration_25antibody.png"),
                p1, device = "png", width = 5, height = 5, units = "in")

############

pbmc2[["ADT"]]@scale.data <- pbmc2[["ADT"]]@scale.data[which(rownames(pbmc2[["ADT"]]@scale.data) %in% antibody_vec2),]

print(dim(pbmc2[["ADT"]]@scale.data))
pbmc2 <- Seurat::RunPCA(
  pbmc2, 
  reduction.name = 'apca',
  features = rownames(pbmc2[["ADT"]]@scale.data),
  verbose = F
)

# compute new ADT umap
print("Computing ADT UMAP")
set.seed(10)
pbmc2 <- Seurat::RunUMAP(
  pbmc2, 
  reduction = 'apca', 
  dims = 1:25, assay = 'ADT', 
  reduction.name = 'adt.umap', 
  reduction.key = 'adtUMAP_'
)

# compute new multimodal neighbor
print("Computing WNN")
pbmc2 <- Seurat::FindMultiModalNeighbors(
  pbmc2, 
  reduction.list = list("pca", "apca"), 
  dims.list = list(1:rank_1, 1:9), 
  modality.weight.name = "RNA.weight"
)

# compute new WNN UMAP
print("Computing WNN UMAP")
set.seed(10)
pbmc2 <- Seurat::RunUMAP(
  pbmc2, 
  nn.name = "weighted.nn", 
  reduction.name = "wnn.umap", 
  reduction.key = "wnnUMAP_"
)

# compute consensus PCA's UMAP
print("Computing Consensus PCA")
n <- ncol(pbmc2)
svd_1 <- svd(pbmc2[["pca"]]@cell.embeddings)
embedding_1 <- multiomicCCA:::.mult_mat_vec(svd_1$u, svd_1$d/svd_1$d[1]*sqrt(n))
svd_2 <- svd(pbmc2[["apca"]]@cell.embeddings)
embedding_2 <- multiomicCCA:::.mult_mat_vec(svd_2$u, svd_2$d/svd_2$d[1]*sqrt(n))
embedding_all <- cbind(embedding_1, embedding_2)
rownames(embedding_all) <- colnames(pbmc2)
pca_res <- stats::prcomp(embedding_all, center = TRUE, scale. = FALSE)
consensus_mat <- pca_res$x[,1:min(rank_1, rank_2)] 

set.seed(10)
consensus_umap <- Seurat::RunUMAP(consensus_mat, 
                                  reduction.key = "umapConsensusPCA_")
pbmc2[["consensus.umap"]] <- Seurat::CreateDimReducObject(consensus_umap@cell.embeddings,
                                                          assay = "ADT")

# rotate all the UMAPs
print("Rotating UMAPs")
anchor_name <- "rna.umap"
other_names <- c("adt.umap", "wnn.umap", "consensus.umap")

for(umap_name in other_names){
  print(umap_name)
  u_mat1 <- pbmc2[[anchor_name]]@cell.embeddings
  u_mat2 <- pbmc2[[umap_name]]@cell.embeddings
  tmp <- svd(t(u_mat1) %*% u_mat2)
  rotation_mat <- tmp$u %*% t(tmp$v)
  tmp <- u_mat2 %*% t(rotation_mat)
  rownames(tmp) <- rownames(pbmc2@meta.data)
  colnames(tmp) <- colnames(pbmc2[[umap_name]]@cell.embeddings)
  pbmc2[[umap_name]]@cell.embeddings <- tmp
}

# make all the plots
print("Plotting")
reduction_vec <- c(other_names)
group_vec <- c("celltype.l2")
main_vec <- c("(ADT)", "(WNN)", "(Consensus PCA)")
file_vec <- c("adt", "wnn", "consensuspca")

for(i in 1:length(reduction_vec)){
  for(j in 1:length(group_vec)){
    plot1 <- Seurat::DimPlot(pbmc2, reduction = reduction_vec[i],
                             group.by = group_vec[j], label = TRUE,
                             repel = TRUE, label.size = 2.5)
    plot1 <- plot1 + ggplot2::ggtitle(paste0("Human PBMC (Cite-seq):\n", main_vec[i]))
    plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
    ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14k/Writeup14k_citeseq_pbmc224_dcca_25antibody_", file_vec[i], ".png"),
                    plot1, device = "png", width = 6, height = 5, units = "in")
  }
}

