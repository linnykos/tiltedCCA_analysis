rm(list=ls()); set.seed(10)

library(Seurat); library(Signac); library(multiomicCCA)

load("../../../../out/Writeup14/Writeup14_10x_mouseembryo_preprocess.RData")
date_of_run <- Sys.time(); session_info <- sessionInfo()

Seurat::DefaultAssay(mbrain) <- "SCT"
mat_1 <- Matrix::t(mbrain[["SCT"]]@data[Seurat::VariableFeatures(object = mbrain),])
Seurat::DefaultAssay(mbrain) <- "ATAC"
mat_2 <- Matrix::t(mbrain[["ATAC"]]@data[Seurat::VariableFeatures(object = mbrain),])

head(rownames(mat_1)); head(colnames(mat_1))
head(rownames(mat_2)); head(colnames(mat_2))
dim(mat_1); dim(mat_2)
metadata <- mbrain@meta.data

cell_idx <- which(metadata$label_Savercat %in% c("Oligodendrocyte", "Radial glia", 
                                                 "Forebrain GABAergic", "Neuroblast", 
                                                 "Glioblast", "Cortical or hippocampal glutamatergic"))

######################

mat_1b <- mat_1[cell_idx,]
sd_vec <- sparseMatrixStats::colSds(mat_1b)
if(any(sd_vec <= 1e-6)){
  mat_1b <- mat_1b[,-which(sd_vec <= 1e-6)]
}

mat_2b <- mat_2[cell_idx,]
sd_vec <- sparseMatrixStats::colSds(mat_2b)
if(any(sd_vec <= 1e-6)){
  mat_2b <- mat_2b[,-which(sd_vec <= 1e-6)]
}

set.seed(10)
rank_1 <- 30; rank_2 <- 50; nn <- 15
dcca_res <- multiomicCCA::dcca_factor(mat_1b, mat_2b, dims_1 = 1:rank_1, dims_2 = 2:rank_2,
                                      center_1 = T, center_2 = T,
                                      meta_clustering = NA, num_neigh = nn, 
                                      fix_distinct_perc = F, verbose = T) 
print("Finished Tilted-CCA")

#######################3

nn = 30
radius_quantile = 0.5
verbose = T
frnn_approx = 0
embedding <- multiomicCCA:::.prepare_embeddings(dcca_res, 
                                                data_1 = T, 
                                                data_2 = F, 
                                                add_noise = F, 
                                                center = T, 
                                                renormalize = F)
n <- nrow(embedding[[1]])

vec_print <- c("common", "distinct", "everything")
vec_rad <- sapply(1:3, function(i){
  if(verbose) print(paste0(Sys.time(),": cLISI: Computing radius -- ", vec_print[i]))
  multiomicCCA:::.compute_radius(embedding[[i]], nn, radius_quantile)
})
vec_rad_org <- vec_rad
names(vec_rad_org) <- c("common", "distinct", "everything")
vec_rad[1:2] <- max(vec_rad[1:2])

# construct frnn
list_g <- lapply(1:3, function(i){
  if(verbose) print(paste0(Sys.time(),": cLISI: Construct graph -- ", vec_print[i]))
  multiomicCCA:::.construct_frnn(embedding[[i]], 
                                 radius = vec_rad[i], 
                                 nn = nn, 
                                 frnn_approx = frnn_approx, 
                                 verbose = verbose)
}) 

for(i in 1:3){
  list_g[[i]] <- multiomicCCA:::.nnlist_to_matrix(list_g[[i]])
  if(symmetrize){
    list_g[[i]] <- multiomicCCA:::.symmetrize_sparse(list_g[[i]], set_ones = F)
  }
  
  # convert back to list form if needed
  if(bool_matrix){
    if(length(rownames(obj$common_score)) != 0){
      rownames(list_g[[i]]) <- rownames(obj$common_score)
      colnames(list_g[[i]]) <- rownames(obj$common_score)
    }
  } else {
    # [[note to self: add a test to make sure this conversion is bijective]]
    list_g[[i]] <- multiomicCCA:::.matrix_to_nnlist(list_g[[i]])
  }
}
