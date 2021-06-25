rm(list=ls()); set.seed(10)

library(Seurat); library(Signac); library(EnsDb.Hsapiens.v86); 
library(BSgenome.Hsapiens.UCSC.hg38); library(GenomeInfoDb)
library(dplyr); library(ggplot2); library(multiomicCCA)

load("../../../../out/Writeup14/Writeup14_10x_mouseembryo_preprocess.RData")
source("Writeup14_laplacian_function")
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
dcca_res <- multiomicCCA::dcca_factor(mat_1[cell_idx,], mat_2[cell_idx,], 
                                      dims_1 = 1:rank_1, dims_2 = 2:rank_2,
                                      center_1 = T, center_2 = T,
                                      meta_clustering = NA, num_neigh = 15, 
                                      apply_shrinkage = F, fix_distinct_perc = F, 
                                      verbose = T) 
dcca_res2 <- dcca_res
class(dcca_res2) <- "dcca_decomp"

dcca_decomp <- multiomicCCA::dcca_decomposition(dcca_res, verbose = T)

#####################
mat_1 <- mat_1[cell_idx,]
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
                     colname_vec = paste0("clap_", 1:200))
d_eig <- compute_lap(rna_frnn$d_g, k_max = 200, rowname_vec = colnames(mbrain2), 
                     colname_vec = paste0("dlap_", 1:200))
e_eig <- compute_lap(rna_frnn$e_g, k_max = 200, rowname_vec = colnames(mbrain2), 
                     colname_vec = paste0("elap_", 1:200))

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
df_tmp <- cbind(mat_1[cell_idx,gene_idx], e_eig)
colnames(df_tmp) <- c("gene", paste0("basis", 1:ncol(e_eig)))
df_tmp <- as.data.frame(df_tmp)
lm_res1 <- stats::lm(gene ~ 1, data = df_tmp)
lm_res2 <- stats::lm(gene ~ ., data = df_tmp)
sqrt(sum(mat_1[cell_idx,gene_idx]^2))
sqrt(sum(lm_res2$fitted.values^2))
stats::var(mat_1[cell_idx,gene_idx])
stats::var(lm_res2$fitted.values)
summary(lm_res1)
anova(lm_res1, lm_res2)
stats::cor(lm_res$fitted.values, mat_1[cell_idx,gene_idx])^2

df_tmp <- cbind(dcca_decomp$common_mat_1[,gene_idx], c_eig)
colnames(df_tmp) <- c("gene", paste0("basis", 1:ncol(c_eig)))
df_tmp <- as.data.frame(df_tmp)
lm_res1 <- stats::lm(gene ~ 1, data = df_tmp)
lm_res2 <- stats::lm(gene ~ ., data = df_tmp)
sqrt(sum(dcca_decomp$common_mat_1[,gene_idx]^2))
sqrt(sum(lm_res2$fitted.values^2))
stats::var(dcca_decomp$common_mat_1[,gene_idx])
stats::var(lm_res2$fitted.values)
summary(lm_res1)
anova(lm_res1, lm_res2)

df_tmp <- cbind(dcca_decomp$distinct_mat_1[,gene_idx], d_eig)
colnames(df_tmp) <- c("gene", paste0("basis", 1:ncol(d_eig)))
df_tmp <- as.data.frame(df_tmp)
lm_res1 <- stats::lm(gene ~ 1, data = df_tmp)
lm_res2 <- stats::lm(gene ~ ., data = df_tmp)
sqrt(sum(dcca_decomp$distinct_mat_1[,gene_idx]^2))
sqrt(sum(lm_res2$fitted.values^2))
stats::var(dcca_decomp$distinct_mat_1[,gene_idx])
stats::var(lm_res2$fitted.values)
summary(lm_res1)
anova(lm_res1, lm_res2)

##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################


set.seed(10)
rna_frnn <- multiomicCCA::construct_frnn(dcca_res$common_score, dcca_res$distinct_score_2, 
                                         svd_e = dcca_res$svd_2, cell_subidx = 1:nrow(dcca_res$common_score), 
                                         nn = 15, bool_matrix = T, verbose = T)

#compute all the degree vectors
c_eig <- compute_lap(rna_frnn$c_g)
d_eig <- compute_lap(rna_frnn$d_g)
e_eig <- compute_lap(rna_frnn$e_g)

#########


rownames(c_eig) <- colnames(mbrain2)
colnames(c_eig) <- paste0("clap_", 1:ncol(c_eig))

set.seed(10)
zz1 <- multiomicCCA::plot_embeddings(dcca_res2, membership_vec, data_1 = F, data_2 = T, 
                                     add_noise = F, pca = F, only_embedding = T, verbose = T)

mbrain2[["common"]] <- Seurat::CreateDimReducObject(embedding = zz1[[1]], key = "UMAP", assay = "RNA")
mbrain2[["clap"]] <- Seurat::CreateDimReducObject(embedding = c_eig, key = "clap", assay = "RNA")

plot1 <- Seurat::FeaturePlot(mbrain2, features = paste0("clap_", 1:16), reduction = "common")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14/Writeup14_10x_mouseembryo_dcca_atac_common_basis.png"),
                plot1, device = "png", width = 16, height = 12, units = "in")

####################

rownames(d_eig) <- colnames(mbrain2)
colnames(d_eig) <- paste0("dlap_", 1:ncol(d_eig))

mbrain2[["distinct"]] <- Seurat::CreateDimReducObject(embedding = zz1[[2]], 
                                                      key = "UMAP", assay = "RNA")
mbrain2[["dlap"]] <- Seurat::CreateDimReducObject(embedding = d_eig, key = "dlap", 
                                                  assay = "RNA")

plot1 <- Seurat::FeaturePlot(mbrain2, features = paste0("dlap_", 1:16), 
                             reduction = "distinct")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14/Writeup14_10x_mouseembryo_dcca_atac_distinct_basis.png"),
                plot1, device = "png", width = 16, height = 12, units = "in")

####################

rownames(e_eig) <- colnames(mbrain2)
colnames(e_eig) <- paste0("elap_", 1:ncol(e_eig))

mbrain2[["everything"]] <- Seurat::CreateDimReducObject(embedding = zz1[[3]], 
                                                        key = "UMAP", assay = "RNA")
mbrain2[["elap"]] <- Seurat::CreateDimReducObject(embedding = e_eig, key = "elap", 
                                                  assay = "RNA")
plot1 <- Seurat::FeaturePlot(mbrain2, features = paste0("elap_", 1:16), reduction = "everything")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14/Writeup14_10x_mouseembryo_dcca_atac_everything_basis.png"),
                plot1, device = "png", width = 16, height = 12, units = "in")



