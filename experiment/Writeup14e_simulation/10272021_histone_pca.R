rm(list=ls())
library(Seurat)
library(Signac)

# meant to be run on desktop
load("../../out/simulation/Writeup14e_simulation/10272021_H3K27me3_original.RData")

cell_list <- list(Non_neurons = c("Astro_Myoc", "Astro_Nnat", "OPC", "Microglia", 
                                  "Oligo_MFOL", "Oligo_MOL", "Endothelial", "Ependymal"),
                  Inhibitory_neurons = c("Pvalb", "Sst", "CGE"),
                  Cortical_neurons = c("CT", "NP", "L6", "L5", "L4", "L23", "PT"),
                  Hippocampal_neurons = c("DG", "Subiculum", "CA1", "CA23"))
color_name_vec <- c("Orange", "Purple", "Blue", "Green")
celltype_custom <- celltype[,1]
celltype_mapping <- color_df$celltype
for(i in 1:length(cell_list)){
  for(j in 1:length(cell_list[[i]])){
    celltype_custom[which(celltype_custom == cell_list[[i]][j])] <- paste0(color_name_vec[i], " ", j)
    celltype_mapping[which(celltype_mapping == cell_list[[i]][j])] <- paste0(color_name_vec[i], " ", j)
  }
}
color_df$custom <- celltype_mapping

mat_1b <- mat_1; colnames(mat_1b) <- paste0("gene-", 1:ncol(mat_1b))
seurat_obj <- Seurat::CreateSeuratObject(counts = t(mat_1b))
mat_2b <- mat_2; colnames(mat_2b) <- paste0("dna-", 1:ncol(mat_2b))
seurat_obj[["DNA"]] <- Seurat::CreateAssayObject(counts = t(mat_2b))
seurat_obj[["pca"]] <- Seurat::CreateDimReducObject(mat_1, assay = "RNA")
seurat_obj[["lsi"]] <- Seurat::CreateDimReducObject(mat_2, assay = "DNA")
seurat_obj[["celltype"]] <- celltype[,1]
seurat_obj[["celltype_custom"]] <- celltype_custom
celltype_l2 <- celltype[,1]
for(i in 1:length(cell_list)){
  celltype_l2[celltype_l2 %in% cell_list[[i]]] <- names(cell_list)[i]
}
seurat_obj[["celltype_l2"]] <- celltype_l2

uniq_celltypes <- sort(unique(seurat_obj@meta.data$celltype_custom))
color_vec <- sapply(uniq_celltypes, function(i){
  color_df[which(color_df$custom == i),"color"]
})

###############################3

# let's start tinkering
pull_factor <- 0.2
mat_2_new <- mat_2
seurat_obj2 <- seurat_obj

celltype_vec <- c("CT", "NP", "L6", "L5", "L4", "L23", "PT")
offset_vec <- c(-.1, .1 , -1.8, 2, 1.1, -4, 2)
for(kk in 1:length(celltype_vec)){
  celltype <- celltype_vec[kk]
  idx <- which(seurat_obj2$celltype == celltype)
  sd_val <- median(matrixStats::colSds(mat_2_new[idx,,drop = F]))/2
  mean_vec <- Matrix::colMeans(mat_2_new[idx,,drop = F])
  for(i in idx){
    mat_2_new[i,] <- pull_factor*mean_vec + (1-pull_factor)*mat_2_new[i,]
    mat_2_new[i,] <- mat_2_new[i,] + rnorm(ncol(mat_2_new), mean = 0, sd = sd_val)
  }
  mat_2_new[idx,1:5] <- mat_2_new[idx,1:5] + offset_vec[kk]
}

idx <- which(seurat_obj2$celltype_l2 != "Cortical_neurons")
mean_vec <- Matrix::colMeans(mat_1_new[idx,,drop = F])
sd_val <- median(matrixStats::colSds(mat_1_new[idx,,drop = F]))/4
for(i in idx){
  mat_1_new[i,] <- pull_factor*mean_vec + (1-pull_factor)*mat_1_new[i,]
  mat_1_new[i,] <- mat_1_new[i,] + rnorm(ncol(mat_1_new), mean = 0, sd = sd_val)
}

mat_2_new <- scale(mat_2_new, center = T, scale = F)
svd_res <- svd(mat_2_new)
svd_res$d[1] <- 27
svd_res$d[2:10] <- exp(seq(log(20), log(10), length.out = 9))
svd_res$d[11:50] <- exp(seq(log(10), log(1), length.out = 40))
mat_2_new <- multiomicCCA:::.mult_mat_vec(svd_res$u, svd_res$d)
colnames(mat_2_new) <- colnames(mat_2)[1:ncol(mat_2_new)]
rownames(mat_2_new) <- rownames(mat_2)
seurat_obj2[["lsi"]] <- Seurat::CreateDimReducObject(mat_2_new, assay = "DNA")

########################################


mat_1_new <- mat_1
rowsum_vec <- apply(mat_1_new, 1, multiomicCCA:::.l2norm)
mat_1_new <- multiomicCCA:::.mult_vec_mat(1/rowsum_vec, mat_1_new)
pull_factor <- 0.7

set.seed(10)
idx <- which(seurat_obj2$celltype_l2 == "Cortical_neurons")
mean_vec <- Matrix::colMeans(mat_1_new[idx,,drop = F])
sd_val <- median(matrixStats::colSds(mat_1_new[idx,,drop = F]))/6
for(i in idx){
  mat_1_new[i,] <- pull_factor*mean_vec + (1-pull_factor)*mat_1_new[i,]
  mat_1_new[i,] <- mat_1_new[i,] + rnorm(ncol(mat_1_new), mean = 0, sd = sd_val)
}

mat_1_new <- scale(mat_1_new, center = T, scale = F)
svd_res <- svd(mat_1_new)
svd_res$d[1:10] <- seq(27, 15, length.out = 10)[1:10]
mat_1_new <- multiomicCCA:::.mult_mat_vec(svd_res$u, svd_res$d)
colnames(mat_1_new) <- colnames(mat_2)[1:ncol(mat_1_new)]
rownames(mat_1_new) <- rownames(mat_2)
seurat_obj2[["pca"]] <- Seurat::CreateDimReducObject(mat_1_new, assay = "RNA")

#############################

mat1 <- seurat_obj2[["pca"]]@cell.embeddings
mat2 <- seurat_obj2[["lsi"]]@cell.embeddings
vec1 <- seurat_obj2[["pca"]]@cell.embeddings[,1]
vec2 <- seurat_obj2[["lsi"]]@cell.embeddings[,1]

pca_res <- stats::prcomp(cbind(vec1, vec2))

xlim <- c(-.8, 0.8) # quantile(c(vec1, vec2), probs = c(0.05, 0.95))

col_vec <- rep(NA, length(vec1))
uniq_celltypes <- sort(unique(seurat_obj2@meta.data$celltype))
for(celltype in uniq_celltypes){
  col_vec[which(seurat_obj2@meta.data$celltype == celltype)] <- color_df[which(color_df$celltype == celltype),"color"]
}
rgb(0.5, 0.5, 0.5, 0.1)
col_vec_trans <- sapply(col_vec, function(x){
  paste0(x, "1A")
})

png("../../out/simulation/Writeup14e_simulation/consensus_leadingpc.png",
    height = 1500, width = 1500, units = "px", res = 300)
par(mar = c(5,5,0.5,0.5))
plot(NA, asp = T, 
     xlim = xlim, ylim = xlim, bty = "n",
     xlab = expression(Leading ~ PC ~ of ~ Modality ~ 1: ~ sigma^(1)==20),
     ylab = expression(Leading ~ PC ~ of ~ Modality ~ 2: ~ sigma^(2)==20))
lines(c(-1e5,1e5), rep(0,2),lty = 2)
lines(rep(0,2), c(-1e5,1e5), lty = 2)
points(vec1, vec2, col = col_vec, pch = 16,)
graphics.off()

n <- length(vec1)
set.seed(10)
jitter_vec <- runif(n, min = -.1, max = .1)
jitter_2d <- jitter_vec %*% t(pca_res$rotation[,2])
projection_vec <- cbind(vec1, vec2) %*% pca_res$rotation[,1] %*% t(pca_res$rotation[,1])
projection_vec <- projection_vec + jitter_2d
png("../../out/simulation/Writeup14e_simulation/consensus_leadingpc_projected.png",
    height = 1500, width = 1500, units = "px", res = 300)
par(mar = c(5,5,0.5,0.5))
plot(NA, asp = T, 
     xlim = xlim, ylim = xlim, bty = "n",
     xlab = expression(Leading ~ PC ~ of ~ Modality ~ 1: ~ sigma^(1)==20),
     ylab = expression(Leading ~ PC ~ of ~ Modality ~ 2: ~ sigma^(2)==20))
lines(c(-1e5,1e5), rep(0,2),lty = 2)
lines(rep(0,2), c(-1e5,1e5), lty = 2)
points(vec1, vec2, col = col_vec_trans, pch = 16)
points(projection_vec, pch = 16, col = "white", cex = 2)
points(projection_vec, pch = 16, col = col_vec)
graphics.off()


png("../../out/simulation/Writeup14e_simulation/example_pcs.png",
    height = 1200, width = 600, units = "px", res = 300)
par(mar = c(0.5, 4, 0.5, 0.5))
plot(NA, ylim = c(-1.2,0.8), xlim = c(-.5, 1.5), 
     ylab = "Leading PCs", xlab = "", bty = "n", xaxt = "n")
points(y = vec1, x = rep(0, length(vec1)) + runif(length(vec1), min = -.3, max = .3),
       pch = 16, col = col_vec, cex = 0.25)
points(y = vec2, x = rep(1, length(vec2)) + runif(length(vec1), min = -.3, max = .3), 
       pch = 16, col = col_vec, cex = 0.25)
graphics.off()
