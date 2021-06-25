rm(list=ls()); set.seed(10)

library(Seurat); library(Signac); library(EnsDb.Hsapiens.v86); 
library(BSgenome.Hsapiens.UCSC.hg38); library(GenomeInfoDb)
library(dplyr); library(ggplot2); library(multiomicCCA)

load("../../../../out/Writeup14/Writeup14_10x_mouseembryo_preprocess.RData")
source("Writeup14_laplacian_function.R")
date_of_run <- Sys.time(); session_info <- sessionInfo()

# Seurat::DefaultAssay(mbrain) <- "ATAC"
# set.seed(10)
# mbrain <- Seurat::ScaleData(mbrain, features = Seurat::VariableFeatures(object = mbrain))

mat_1 <- t(mbrain[["SCT"]]@scale.data)
Seurat::DefaultAssay(mbrain) <- "ATAC"
mat_2 <- Matrix::t(mbrain[["ATAC"]]@data[Seurat::VariableFeatures(object = mbrain),])

head(rownames(mat_1)); head(colnames(mat_1))
head(rownames(mat_2)); head(colnames(mat_2))
dim(mat_1); dim(mat_2)
metadata <- mbrain@meta.data

cell_idx <- which(metadata$label_Savercat %in% c("Oligodendrocyte", "Hindbrain glycinergic", "Midbrain glutamatergic",
                                                 "Forebrain glutamatergic", "Radial glia", "Forebrain GABAergic", 
                                                 "Neuroblast", "Cajal-Retzius", "Mixed region GABAergic", 
                                                 "Glioblast", "Cortical or hippocampal glutamatergic"))
set.seed(10)
rank_1 <- 30; rank_2 <- 31
mat_1 <- mat_1[cell_idx,]; mat_2 <- mat_2[cell_idx,]
dcca_res <- multiomicCCA::dcca_factor(mat_1, mat_2, dims_1 = 1:rank_1, dims_2 = 2:rank_2,
                                      center_1 = T, center_2 = T,
                                      meta_clustering = NA, num_neigh = 15, 
                                      apply_shrinkage = F, fix_distinct_perc = F, 
                                      verbose = T) 
dcca_res2 <- dcca_res
class(dcca_res2) <- "dcca_decomp"

dcca_decomp <- multiomicCCA::dcca_decomposition(dcca_res, verbose = T)
mat_1_denoised <- dcca_decomp$common_mat_1 + dcca_decomp$distinct_mat_1

#####################

mbrain2 <- Seurat::CreateSeuratObject(counts = t(mat_1))
mbrain2[["label_Savercat"]] <- metadata$label_Savercat[cell_idx]
membership_vec <- as.factor(metadata$label_Savercat[cell_idx])

set.seed(10)
rna_frnn <- multiomicCCA::construct_frnn(dcca_res$common_score, dcca_res$distinct_score_1, 
                           svd_e = dcca_res$svd_1, cell_subidx = 1:nrow(dcca_res$common_score), 
                           nn = 15, bool_matrix = T, include_diag = F, verbose = T)

set.seed(10)
zz1 <- multiomicCCA::plot_embeddings(dcca_res2, membership_vec, data_1 = T, data_2 = F, 
                                     add_noise = F, pca = F, only_embedding = T, verbose = T)

#compute all the degree vectors
c_eig <- compute_lap(rna_frnn$c_g, k_max = 200, rowname_vec = colnames(mbrain2), 
                     colname_vec = paste0("clap_", 1:199))
d_eig <- compute_lap(rna_frnn$d_g, k_max = 200, rowname_vec = colnames(mbrain2), 
                     colname_vec = paste0("dlap_", 1:199))
e_eig <- compute_lap(rna_frnn$e_g, k_max = 200, rowname_vec = colnames(mbrain2), 
                     colname_vec = paste0("elap_", 1:199))

#########

mbrain2[["common"]] <- Seurat::CreateDimReducObject(embedding = zz1[[1]], key = "UMAP", assay = "RNA")
mbrain2[["clap"]] <- Seurat::CreateDimReducObject(embedding = c_eig, key = "clap", assay = "RNA")

plot1 <- Seurat::FeaturePlot(mbrain2, features = paste0("clap_", 1:16), reduction = "common")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14/Writeup14_10x_mouseembryo_dcca_rna_common_basis.png"),
                plot1, device = "png", width = 16, height = 12, units = "in")

###

mbrain2[["distinct"]] <- Seurat::CreateDimReducObject(embedding = zz1[[2]], key = "UMAP", assay = "RNA")
mbrain2[["dlap"]] <- Seurat::CreateDimReducObject(embedding = d_eig, key = "dlap", assay = "RNA")

plot1 <- Seurat::FeaturePlot(mbrain2, features = paste0("dlap_", 1:16), 
                             reduction = "distinct")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14/Writeup14_10x_mouseembryo_dcca_rna_distinct_basis.png"),
                plot1, device = "png", width = 16, height = 12, units = "in")

###

mbrain2[["everything"]] <- Seurat::CreateDimReducObject(embedding = zz1[[3]], key = "UMAP", assay = "RNA")
mbrain2[["elap"]] <- Seurat::CreateDimReducObject(embedding = e_eig, key = "elap", assay = "RNA")

plot1 <- Seurat::FeaturePlot(mbrain2, features = paste0("elap_", 1:16), reduction = "everything")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14/Writeup14_10x_mouseembryo_dcca_rna_everything_basis.png"),
                plot1, device = "png", width = 16, height = 12, units = "in")

########

gene_idx <- which(colnames(mat_1) == "Neurod2")
c_res <- compute_smooth_signal(dcca_decomp$common_mat_1[,gene_idx], c_eig)
d_res <- compute_smooth_signal(dcca_decomp$distinct_mat_1[,gene_idx], d_eig)
e_res <- compute_smooth_signal(mat_1_denoised[,gene_idx], e_eig)
var(e_res$pred_vec); var(c_res$pred_vec); var(d_res$pred_vec)

mbrain3 <- mbrain2
create_plot(mbrain3, var_name = "Neurod2", e_vec = mat_1_denoised[,gene_idx],
            c_vec = dcca_decomp$common_mat_1[,gene_idx],
            d_vec = dcca_decomp$distinct_mat_1[,gene_idx],
            e_res = e_res, c_res = c_res, d_res = d_res, 
            filename = paste0("../../../../out/figures/Writeup14/Writeup14_10x_mouseembryo_dcca_rna_Neurod2.png"))




##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
mbrain2 <- Seurat::CreateSeuratObject(counts = t(mat_1))
mbrain2[["label_Savercat"]] <- metadata$label_Savercat[cell_idx]
membership_vec <- as.factor(metadata$label_Savercat[cell_idx])

set.seed(10)
atac_frnn <- multiomicCCA::construct_frnn(dcca_res$common_score, dcca_res$distinct_score_2, 
                                         svd_e = dcca_res$svd_2, cell_subidx = 1:nrow(dcca_res$common_score), 
                                         nn = 15, bool_matrix = T, include_diag = F, verbose = T)

set.seed(10)
zz2 <- multiomicCCA::plot_embeddings(dcca_res2, membership_vec, data_1 = F, data_2 = T, 
                                     add_noise = F, pca = F, only_embedding = T, verbose = T)

#compute all the degree vectors
c_eig <- compute_lap(atac_frnn$c_g, k_max = 200, rowname_vec = colnames(mbrain2), 
                     colname_vec = paste0("clap_", 1:199))
d_eig <- compute_lap(atac_frnn$d_g, k_max = 200, rowname_vec = colnames(mbrain2), 
                     colname_vec = paste0("dlap_", 1:199))
e_eig <- compute_lap(atac_frnn$e_g, k_max = 200, rowname_vec = colnames(mbrain2), 
                     colname_vec = paste0("elap_", 1:199))


#########

mbrain2[["common"]] <- Seurat::CreateDimReducObject(embedding = zz2[[1]], key = "UMAP", assay = "RNA")
mbrain2[["clap"]] <- Seurat::CreateDimReducObject(embedding = c_eig, key = "clap", assay = "RNA")

plot1 <- Seurat::FeaturePlot(mbrain2, features = paste0("clap_", 1:16), reduction = "common")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14/Writeup14_10x_mouseembryo_dcca_atac_common_basis.png"),
                plot1, device = "png", width = 16, height = 12, units = "in")

####################

mbrain2[["distinct"]] <- Seurat::CreateDimReducObject(embedding = zz2[[2]], key = "UMAP", assay = "RNA")
mbrain2[["dlap"]] <- Seurat::CreateDimReducObject(embedding = d_eig, key = "dlap",  assay = "RNA")

plot1 <- Seurat::FeaturePlot(mbrain2, features = paste0("dlap_", 1:16), reduction = "distinct")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14/Writeup14_10x_mouseembryo_dcca_atac_distinct_basis.png"),
                plot1, device = "png", width = 16, height = 12, units = "in")

####################

mbrain2[["everything"]] <- Seurat::CreateDimReducObject(embedding = zz2[[3]], key = "UMAP", assay = "RNA")
mbrain2[["elap"]] <- Seurat::CreateDimReducObject(embedding = e_eig, key = "elap",  assay = "RNA")

plot1 <- Seurat::FeaturePlot(mbrain2, features = paste0("elap_", 1:16), reduction = "everything")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14/Writeup14_10x_mouseembryo_dcca_atac_everything_basis.png"),
                plot1, device = "png", width = 16, height = 12, units = "in")



