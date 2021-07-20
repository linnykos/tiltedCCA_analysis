rm(list=ls()); set.seed(10)

library(Seurat); library(Signac); library(EnsDb.Hsapiens.v86); 
library(BSgenome.Hsapiens.UCSC.hg38); library(GenomeInfoDb)
library(dplyr); library(ggplot2); library(multiomicCCA)

load("../../../../out/Writeup10/Writeup10_10x_pbmc_preprocess4.RData")
date_of_run <- Sys.time(); session_info <- sessionInfo()

tab <- table(pbmc@meta.data$predicted.id)
cell_idx <- which(pbmc@meta.data$predicted.id %in% names(tab)[tab > 5])

mat_1 <- t(pbmc[["SCT"]]@scale.data[,cell_idx])
mat_2 <- t(pbmc[["ATAC"]]@scale.data[,cell_idx])

head(rownames(mat_1)); head(colnames(mat_1))
head(rownames(mat_2)); head(colnames(mat_2))
dim(mat_1); dim(mat_2)
metadata <- pbmc@meta.data[cell_idx,]
table(metadata$predicted.id)
rm(list = "pbmc"); gc(T)

# set.seed(10)
# rank_1 <- 30; rank_2 <- 50; nn <- 30
# dcca_res <- multiomicCCA::dcca_factor(mat_1, mat_2, dims_1 = 1:rank_1, dims_2 = 2:rank_2,
#                                       meta_clustering = NA, num_neigh = nn, 
#                                       apply_shrinkage = F, fix_distinct_perc = F, 
#                                       verbose = T) 
load("../../../../out/Writeup14b/Writeup14b_10x_pbmc_dcca_laplacian_variablecalculations.RData"); nn <- 30

dcca_decomp <- multiomicCCA::dcca_decomposition(dcca_res, verbose = T)
mat_1_denoised <- dcca_decomp$common_mat_1 + dcca_decomp$distinct_mat_1
mat_2_denoised <- dcca_decomp$common_mat_2 + dcca_decomp$distinct_mat_2

rm(list = c("mat_1", "mat_2")); gc(T)

#########################

pbmc2 <- Seurat::CreateSeuratObject(counts = t(mat_1_denoised))
pbmc2[["celltype"]] <- metadata$predicted.id
membership_vec <- as.factor(metadata$predicted.id)

title_vec <- c("Common view", "Distinct view", "Everything view")
main_vec <- c("common", "distinct", "everything")

set.seed(10)
rna_frnn <- multiomicCCA::construct_frnn(dcca_res, nn = nn, membership_vec = membership_vec,
                                         data_1 = T, data_2 = F,
                                         radius_quantile = 0.5, symmetrize = F, 
                                         bool_matrix = T, verbose = T)

#compute all rna basis vectors
k_max <- 200
c_eig <- multiomicCCA::compute_laplacian(rna_frnn$c_g, k_max = k_max, rowname_vec = colnames(pbmc2), 
                                         colname_vec = paste0("clap_", 1:k_max))
d_eig <- multiomicCCA::compute_laplacian(rna_frnn$d_g, k_max = k_max, rowname_vec = colnames(pbmc2), 
                                         colname_vec = paste0("dlap_", 1:k_max))
e_eig <- multiomicCCA::compute_laplacian(rna_frnn$e_g, k_max = k_max, rowname_vec = colnames(pbmc2), 
                                         colname_vec = paste0("elap_", 1:k_max))

pbmc2[["clap"]] <- Seurat::CreateDimReducObject(embedding = c_eig, key = "clap", assay = "RNA")
pbmc2[["dlap"]] <- Seurat::CreateDimReducObject(embedding = d_eig, key = "dlap", assay = "RNA")
pbmc2[["elap"]] <- Seurat::CreateDimReducObject(embedding = e_eig, key = "elap", assay = "RNA")

set.seed(10)
rna_embeddings <- multiomicCCA::plot_embeddings2(dcca_res, nn = nn, data_1 = T, data_2 = F,
                                                 c_g = rna_frnn$c_g, d_g = rna_frnn$d_g, 
                                                 only_embedding = T, verbose = T, 
                                                 sampling_type = "adaptive_gaussian",
                                                 keep_nn = T)
pbmc2[["common"]] <- Seurat::CreateDimReducObject(embedding = rna_embeddings[[1]], key = "UMAP", assay = "RNA")
pbmc2[["distinct"]] <- Seurat::CreateDimReducObject(embedding = rna_embeddings[[2]], key = "UMAP", assay = "RNA")
pbmc2[["everything"]] <- Seurat::CreateDimReducObject(embedding = rna_embeddings[[3]], key = "UMAP", assay = "RNA")

################

set.seed(10)
atac_frnn <- multiomicCCA::construct_frnn(dcca_res, nn = nn, membership_vec = membership_vec,
                                          data_1 = F, data_2 = T,
                                          radius_quantile = 0.5,
                                          bool_matrix = T, symmetrize = F, verbose = T)

#compute all the degree vectors
k_max <- 200
c_eig2 <- multiomicCCA::compute_laplacian(atac_frnn$c_g, k_max = k_max, rowname_vec = colnames(pbmc2), 
                                          colname_vec = paste0("clap_", 1:k_max))
d_eig2 <- multiomicCCA::compute_laplacian(atac_frnn$d_g, k_max = k_max, rowname_vec = colnames(pbmc2), 
                                          colname_vec = paste0("dlap_", 1:k_max))
e_eig2 <- multiomicCCA::compute_laplacian(atac_frnn$e_g, k_max = k_max, rowname_vec = colnames(pbmc2), 
                                          colname_vec = paste0("elap_", 1:k_max))

pbmc2[["clap2"]] <- Seurat::CreateDimReducObject(embedding = c_eig2, key = "clap", assay = "RNA")
pbmc2[["dlap2"]] <- Seurat::CreateDimReducObject(embedding = d_eig2, key = "dlap", assay = "RNA")
pbmc2[["elap2"]] <- Seurat::CreateDimReducObject(embedding = e_eig2, key = "elap", assay = "RNA")

set.seed(10)
atac_embeddings <- multiomicCCA::plot_embeddings2(dcca_res, nn = nn, data_1 = F, data_2 = T,
                                                  c_g = atac_frnn$c_g, d_g = atac_frnn$d_g, 
                                                  only_embedding = T, verbose = T, 
                                                  sampling_type = "adaptive_gaussian",
                                                  keep_nn = T)
pbmc2[["common2"]] <- Seurat::CreateDimReducObject(embedding = atac_embeddings[[1]], key = "UMAP", assay = "RNA")
pbmc2[["distinct2"]] <- Seurat::CreateDimReducObject(embedding = atac_embeddings[[2]], key = "UMAP", assay = "RNA")
pbmc2[["everything2"]] <- Seurat::CreateDimReducObject(embedding = atac_embeddings[[3]], key = "UMAP", assay = "RNA")

##########################

set.seed(10)
combined_g <- multiomicCCA::combine_frnn(dcca_res, g_1 = rna_frnn$c_g,
                                         g_2 = atac_frnn$c_g, nn = nn)
combined_g <- SeuratObject::as.Graph(combined_g)
set.seed(10)
combined_common_umap <- Seurat::RunUMAP(combined_g, assay = "RNA")@cell.embeddings
pbmc2[["combined"]] <- Seurat::CreateDimReducObject(embedding = combined_common_umap, key = "UMAP", assay = "RNA")

set.seed(10)
both_embeddings <- multiomicCCA::plot_embeddings(dcca_res, membership_vec = membership_vec,
                                                 data_1 = T, data_2 = T, add_noise = F,
                                                 only_embedding = T)
pbmc2[["both"]] <- Seurat::CreateDimReducObject(embedding = both_embeddings$everything, key = "UMAP", assay = "RNA")

##########################################3
# done with all the calculations (for now). now plot

# plot RNA embeddings
for(i in 1:3){
  plot1 <- Seurat::DimPlot(pbmc2, reduction = main_vec[i], 
                           group.by = "celltype", label = TRUE,
                           repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("10x PBMC (RNA)\n", title_vec[i]))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_10x_pbmc_dcca_rna_", main_vec[i], "_umap.png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
}

# plot ATAC embeddings
for(i in 1:3){
  plot1 <- Seurat::DimPlot(pbmc2, reduction = paste0(main_vec[i], "2"), 
                           group.by = "celltype", label = TRUE,
                           repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("10x PBMC (ATAC)\n", title_vec[i])) 
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_10x_pbmc_dcca_atac_", main_vec[i], "_umap.png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
}

# plot the combined view
plot1 <- Seurat::DimPlot(pbmc2, reduction = "combined", 
                         group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle("10x PBMC\nBoth Common") 
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_10x_pbmc_dcca_both_common_umap.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

# plot the everything combined for both
plot1 <- Seurat::DimPlot(pbmc2, reduction = "both", 
                         group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle("10x PBMC\nBoth Everything") 
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_10x_pbmc_dcca_both_everything_umap.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")


#######################################
# now plot the bases -- first RNA
plot1 <- Seurat::FeaturePlot(pbmc2, features = paste0("clap_", 1:16), reduction = "combined")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_10x_pbmc_dcca_rna_common_basis.png"),
                plot1, device = "png", width = 16, height = 12, units = "in")
 
plot1 <- Seurat::FeaturePlot(pbmc2, features = paste0("dlap_", 1:16), reduction = "distinct")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_10x_pbmc_dcca_rna_distinct_basis.png"),
                plot1, device = "png", width = 16, height = 12, units = "in")

plot1 <- Seurat::FeaturePlot(pbmc2, features = paste0("elap_", 1:16), reduction = "both")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_10x_pbmc_dcca_rna_everything_basis.png"),
                plot1, device = "png", width = 16, height = 12, units = "in")

# next ATAC
plot1 <- Seurat::FeaturePlot(pbmc2, features = paste0("clap2_", 1:16), reduction = "combined")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_10x_pbmc_dcca_atac_common_basis.png"),
                plot1, device = "png", width = 16, height = 12, units = "in")

plot1 <- Seurat::FeaturePlot(pbmc2, features = paste0("dlap2_", 1:16), reduction = "distinct2")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_10x_pbmc_dcca_atac_distinct_basis.png"),
                plot1, device = "png", width = 16, height = 12, units = "in")

plot1 <- Seurat::FeaturePlot(pbmc2, features = paste0("elap2_", 1:16), reduction = "both")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_10x_pbmc_dcca_atac_everything_basis.png"),
                plot1, device = "png", width = 16, height = 12, units = "in")

########################


# compute local enrichment
set.seed(10)
rna_local <- multiomicCCA::clisi_information(rna_frnn$c_g, rna_frnn$d_g, rna_frnn$e_g, 
                                             membership_vec = membership_vec)
# rna_local$common_clisi$membership_info
# rna_local$distinct_clisi$membership_info

set.seed(10)
atac_local <- multiomicCCA::clisi_information(atac_frnn$c_g, atac_frnn$d_g, atac_frnn$e_g, 
                                              membership_vec = membership_vec)
# atac_local$common_clisi$membership_info
# atac_local$distinct_clisi$membership_info

tmp <- multiomicCCA::plot_clisi(rna_local, atac_local, main1 = "RNA", main2 = "ATAC")
tmp2 <- cowplot::plot_grid(tmp[[1]], tmp[[2]])
cowplot::save_plot(filename = "../../../../out/figures/Writeup14b/Writeup14b_10x_pbmc_enrichment.png", 
                   tmp2, ncol = 1, nrow = 2, base_height = 1.75, base_asp = 4, device = "png")

#######################



