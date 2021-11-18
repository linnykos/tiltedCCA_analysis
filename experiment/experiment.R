obj = dcca_res2
membership_vec = NA
nn = 30
data_1 = F
data_2 = T
normalization_type = "cosine_itself"
max_subsample_frnn = nrow(obj$common_score)
frnn_approx = 0
radius_quantile = 0.5
bool_matrix = T
center = T
renormalize = F
symmetrize = F
verbose = T

###########

stopifnot(frnn_approx >= 0, frnn_approx <= 1)
if(!all(is.na(membership_vec))){
  stopifnot(length(membership_vec) == nrow(obj$common_score),
            is.factor(membership_vec))
}

embedding <- multiomicCCA:::.prepare_embeddings(obj, data_1 = data_1, data_2 = data_2, 
                                 center = center, 
                                 renormalize = renormalize)
embedding <- multiomicCCA:::.normalize_embeddings(embedding, normalization_type)
n <- nrow(embedding[[1]])

# construct subsamples
if(!all(is.na(membership_vec))){
  cell_subidx <- multiomicCCA:::.construct_celltype_subsample(membership_vec, max_subsample_frnn)
  if(length(cell_subidx) < n) {
    membership_vec <- membership_vec[cell_subidx]
  }
  for(i in 1:3){
    embedding[[i]] <- embedding[[i]][cell_subidx,,drop = F]
  }
}
n <- nrow(embedding[[1]])

######################

vec_print <- c("common", "distinct")
vec_rad <- sapply(1:2, function(i){
  if(verbose) print(paste0(Sys.time(),": cLISI: Computing radius -- ", vec_print[i]))
  multiomicCCA:::.compute_radius(embedding[[i]], nn, radius_quantile)
})
vec_rad_org <- vec_rad
names(vec_rad_org) <- c("common", "distinct")
vec_rad[1:2] <- max(vec_rad[1:2])

#################

list_g <- lapply(1:2, function(i){
  if(verbose) print(paste0(Sys.time(),": cLISI: Construct graph -- ", vec_print[i]))
  multiomicCCA:::.construct_frnn(embedding[[i]], radius = vec_rad[i], nn = nn, 
                  frnn_approx = frnn_approx, 
                  resolve_isolated_nodes = T,
                  radius_quantile = NA,
                  verbose = verbose)
}) 

########################

i = 1
list_g[[i]] <-  multiomicCCA:::.nnlist_to_matrix(list_g[[i]], set_to_one = F)
if(symmetrize){
  list_g[[i]] <- .symmetrize_sparse(list_g[[i]], set_ones = F)
}

# convert back to list form if needed
if(bool_matrix){
  if(length(rownames(obj$common_score)) != 0){
    rownames(list_g[[i]]) <- rownames(obj$common_score)
    colnames(list_g[[i]]) <- rownames(obj$common_score)
  }
} else {
  # [[note to self: add a test to make sure this conversion is bijective]]
  list_g[[i]] <- .matrix_to_nnlist(list_g[[i]])
}

i = 2
list_g[[i]] <-  multiomicCCA:::.nnlist_to_matrix(list_g[[i]], set_to_one = F)
if(symmetrize){
  list_g[[i]] <- .symmetrize_sparse(list_g[[i]], set_ones = F)
}

# convert back to list form if needed
if(bool_matrix){
  if(length(rownames(obj$common_score)) != 0){
    rownames(list_g[[i]]) <- rownames(obj$common_score)
    colnames(list_g[[i]]) <- rownames(obj$common_score)
  }
} else {
  # [[note to self: add a test to make sure this conversion is bijective]]
  list_g[[i]] <- .matrix_to_nnlist(list_g[[i]])
}



