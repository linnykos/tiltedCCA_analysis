rm(list=ls()); set.seed(10)

library(Seurat); library(Signac); library(EnsDb.Hsapiens.v86); 
library(BSgenome.Hsapiens.UCSC.hg38); library(GenomeInfoDb)
library(dplyr); library(ggplot2); library(multiomicCCA)

load("../../../../out/Writeup14/Writeup14_10x_mouseembryo_preprocess.RData")
source("Writeup14_peakcalling_function.R")
date_of_run <- Sys.time(); session_info <- sessionInfo()

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
mat_2_denoised <- dcca_decomp$common_mat_2 + dcca_decomp$distinct_mat_2

############

mbrain2 <- Seurat::CreateSeuratObject(counts = t(mat_1_denoised))
mbrain2[["label_Savercat"]] <- metadata$label_Savercat[cell_idx]
membership_vec <- as.factor(metadata$label_Savercat[cell_idx])

set.seed(10)
rna_frnn <- multiomicCCA::construct_frnn(dcca_res, nn = 15, membership_vec = membership_vec,
                                         data_1 = T, data_2 = F,
                                         bool_matrix = T, include_diag = F, verbose = T)

############3
# starting to debug what's going on:
.distance_to_kernel2 <- function(mat_g){
  stopifnot(inherits(mat_g, "dgCMatrix"))
  
  n <- nrow(mat_g)
  tmp <- lapply(1:n, function(col_idx){
    val1 <- mat_g@p[col_idx]+1
    val2 <- mat_g@p[col_idx+1]+1
    
    if(val1 == val2) return(numeric(0))
    vec <- mat_g@x[(val1+1):min(val2,length(mat_g@x))]
    min_val <- min(vec)
    max_val <- max(vec)
    exp(-(vec - min_val)/max_val)
  })
  print(which(sapply(tmp, function(x){any(is.na(x))})))
  x_val <- unlist(tmp)
  print(length(which(is.na(x_val))))

  mat_g@x <- x_val
  mat_g
}

i <- 1
list_g <- rna_frnn
list_g[[i]] <-  multiomicCCA:::.symmetrize_sparse(list_g[[i]], set_ones = F)
quantile(list_g[[i]]@x)
yy <- list_g[[i]]@x

# convert into kernels
mat_g <- list_g[[i]]
list_g[[i]] <- .distance_to_kernel2(list_g[[i]])
zz <- list_g[[i]]@x
quantile(zz)
length(which(is.na(zz)))
list_g[[i]] <- SeuratObject::as.Graph(list_g[[i]])

set.seed(10)
zz <- Seurat::RunUMAP(list_g[[i]], verbose = T)@cell.embeddings

set.seed(10)
zz1 <- multiomicCCA::plot_embeddings2(c_g = rna_frnn[[1]], d_g = rna_frnn[[2]], e_g = rna_frnn[[3]], 
                                      membership_vec, data_1 = T, data_2 = F, 
                                      only_embedding = T, verbose = T)

#compute all the degree vectors
c_eig <- multiomicCCA::compute_laplacian(rna_frnn$c_g, k_max = 200, rowname_vec = colnames(mbrain2), 
                     colname_vec = paste0("clap_", 1:199))
d_eig <- multiomicCCA::compute_laplacian(rna_frnn$d_g, k_max = 200, rowname_vec = colnames(mbrain2), 
                     colname_vec = paste0("dlap_", 1:199))
e_eig <- multiomicCCA::compute_laplacian(rna_frnn$e_g, k_max = 200, rowname_vec = colnames(mbrain2), 
                     colname_vec = paste0("elap_", 1:199))

############################
title_vec <- c("Common view", "Distinct view", "Everything view")
main_vec <- c("common", "distinct", "everything")

for(i in 1:3){
  mbrain2[["asdf"]] <- Seurat::CreateDimReducObject(embedding = zz, key = "UMAP", assay = "RNA")
  plot1 <- Seurat::DimPlot(mbrain2, reduction = "asdf", group.by = "label_Savercat", label = TRUE,
                           repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle("Mouse Embryo (RNA)\n", title_vec[i]) 
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14/Writeup14_10x_mouseembryo_dcca_rna_", main_vec[i], "_umap2.png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
}
