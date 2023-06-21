generate_B <- function(shrink = 0){
  stopifnot(shrink >= 0, shrink <= 1)
  
  shrink2 <- min(2*shrink, 1)
  val1 <- (1-shrink2)*0.6 + shrink2*0.89
  val2 <- (1-shrink)*0.6 + shrink*0.89
  
  B_mat <- matrix(c(0.9, val1, val2,
                    val1, 0.9, val2,
                    val2, val2, 0.9), 3, 3, byrow = T)
}

simulation_function <- function(shrink = 1){
  n_each <- 100
  B_mat1 <- generate_B(shrink = 0)
  B_mat2 <- generate_B(shrink = shrink)
  K <- ncol(B_mat1)
  membership_vec <- rep(1:3, each = n_each)
  true_cluster <- membership_vec
  n <- length(membership_vec); true_membership_vec <- membership_vec
  svd_u_1 <- tiltedCCA::generate_sbm_orthogonal(B_mat1, membership_vec, centered = T)[,1:2]
  svd_u_2 <- tiltedCCA::generate_sbm_orthogonal(B_mat2, membership_vec, centered = T)[,1:2]
  svd_u_2[,1] <- svd_u_2[,1]/sqrt(sum(svd_u_2[,1]^2))
  svd_u_2[,2] <- svd_u_2[,2]/sqrt(sum(svd_u_2[,2]^2))
  
  p <- 200
  svd_d_1 <- sqrt(n*p)*c(1.5,1); svd_d_2 <- sqrt(n*p)*c(1.5,1)
  svd_v_1 <- tiltedCCA::generate_random_orthogonal(p, 2)
  svd_v_2 <- tiltedCCA::generate_random_orthogonal(p, 2)
  
  mat_1 <- tcrossprod(svd_u_1 %*% diag(svd_d_1), svd_v_1)
  mat_2 <- tcrossprod(svd_u_2 %*% diag(svd_d_2), svd_v_2)
  
  # add some noise
  sd_vec <- c(0.5, 1, 5)
  for(j in 1:p){
    idx <- j %% 3; if(idx == 0) idx <- 3
    mat_1[,j] <- mat_1[,j] + stats::rnorm(nrow(mat_1), mean = 0, sd = sd_vec[idx])
    mat_2[,j] <- mat_2[,j] + stats::rnorm(nrow(mat_2), mean = 0, sd = sd_vec[idx])
  }
  
  mat_1 <- scale(mat_1, center = T, scale = T)
  mat_2 <- scale(mat_2, center = T, scale = T)
  svd_2 <- svd(mat_2)
  
  # number of clusters for mat_2
  mclust_res <- mclust::Mclust(svd_2$u %*% diag(svd_2$d), 
                               G = 1:3,
                               modelNames = "VVV",
                               verbose = F)
  
  clustering_1 <- factor(stats::kmeans(mat_1, centers = 3)$cluster)
  clustering_2 <- factor(stats::kmeans(mat_2, centers = mclust_res$G)$cluster)
  
  rownames(mat_1) <- paste0("n", 1:nrow(mat_1))
  rownames(mat_2) <- paste0("n", 1:nrow(mat_2))
  colnames(mat_1) <- paste0("g", 1:ncol(mat_1))
  colnames(mat_2) <- paste0("p", 1:ncol(mat_2))
  
  list(clustering_1 = clustering_1,
       clustering_2 = clustering_2,
       mat_1 = mat_1,
       mat_2 = mat_2,
       true_cluster = true_cluster)
}

svd_func <- function(mat){
  svd_res <- svd(mat)
  svd_res$u[,1:2] %*% diag(svd_res$d[1:2])
}

compute_alignment <- function(mat_1,
                              mat_2,
                              multiSVD_obj,
                              true_cluster){
  seurat_obj <- Seurat::CreateSeuratObject(counts = t(mat_1))
  seurat_obj[["ADT"]] <- Seurat::CreateAssayObject(counts = t(mat_2))
  seurat_obj$label <- true_cluster
  
  seurat_obj[["RNA"]]@var.features <- rownames(seurat_obj[["RNA"]]@counts)
  seurat_obj[["ADT"]]@var.features <- rownames(seurat_obj[["ADT"]]@counts)
  
  de_list <- tiltedCCA:::differential_expression(seurat_obj = seurat_obj,
                                                 assay = "RNA",
                                                 idents = "label",
                                                 test_use = "wilcox",
                                                 slot = "data",
                                                 verbose = F)
  
  gene_names <- rownames(seurat_obj[["RNA"]]@counts)
  celltype_names <- sort(unique(seurat_obj$label))
  logpval_vec <- sapply(1:length(gene_names), function(k){
    gene <- gene_names[k]
    
    # cycle through all the celltypes
    celltype_vec <- sapply(1:length(de_list$level_vec), function(i){
      idx <- which(de_list$combn_mat == i, arr.ind = T)[,2]
      vec <-  sapply(idx, function(j){
        idx <- which(rownames(de_list$de_list[[j]]) == gene)
        if(length(idx) == 0) return(1)
        de_list$de_list[[j]][idx, "p_val"]
      })
      stats::quantile(vec, probs = 0.75)
    })
    names(celltype_vec) <- celltype_names
    
    max(-log10(celltype_vec))
  })
  names(logpval_vec) <- Seurat::VariableFeatures(seurat_obj)
  
  rsquare_vec <- tiltedCCA:::postprocess_modality_alignment(input_obj = multiSVD_obj,
                                                            bool_use_denoised = T,
                                                            seurat_obj = seurat_obj,
                                                            input_assay = 1,
                                                            seurat_assay = "RNA",
                                                            seurat_slot = "data",
                                                            verbose = F)
  all(names(logpval_vec) == names(rsquare_vec))
  stats::median(rsquare_vec[which(logpval_vec >= 10)])
}