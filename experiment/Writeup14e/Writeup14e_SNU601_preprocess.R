rm(list=ls())
genotype_est <- readRDS("../../../../data/Chiyun_SNU601/SNU601/genotype_est.rds")
genotype_values <- readRDS("../../../../data/Chiyun_SNU601/SNU601/genotype_values_h1_h2_c20.rds")
seg_table_filtered <- readRDS("../../../../data/Chiyun_SNU601/SNU601/seg_table_filtered.rds")
mat <- Matrix::readMM("../../../../data/Chiyun_SNU601/SNU601/matrix.mtx")
features <- read.csv("../../../../data/Chiyun_SNU601/SNU601/features.txt", header = F)
barcodes <- read.csv("../../../../data/Chiyun_SNU601/SNU601/barcodes.txt", header = F)
clone_assign <- readRDS("../../../../data/Chiyun_SNU601/SNU601/cloneAssign.rds")

rownames(mat) <- features[,1]
colnames(mat) <- barcodes[,1]

# remove cells that are not assigned to any clones
cell_names <- names(clone_assign)[which(!is.na(clone_assign))]
if(length(cell_names) > 0){
  mat <- mat[,which(colnames(mat) %in% cell_names)]
  genotype_values <- genotype_values[which(rownames(genotype_values) %in% cell_names),]
  genotype_est <- genotype_est[which(rownames(genotype_est) %in% cell_names),]
  clone_assign <- clone_assign[which(names(clone_assign) %in% cell_names)]
}
mat <- mat[,names(clone_assign)]

######## let's first try using the estimated phi's and theta's

# modality 2 is based on the estimated state. I want both the discrete phi and theta estimated
# [[what does NA mean?]]
# from https://github.com/seasoncloud/Alleloscope/blob/main/R/genotype_ref.R
canonicalPoints <- data.frame(code = c(1:27),
                              rho = c(0.5,0.5, 1,1,1,1.5,1.5,1.5,1.5,2,2,2,2,2,2.5,2.5,2.5,2.5,2.5,2.5,3,3,3,3,3,3,3),
                              theta = c(0,1,0,0.5,1,0,1/3,2/3,1,0,1/4,2/4,3/4,1,0,1/5,2/5,3/5,4/5,1,0,1/6,2/6,3/6,4/6,5/6,1))
tmp <- genotype_est
mat_rho <- tmp
for(i in 1:27){
  mat_rho[tmp == i] <- canonicalPoints[i,"rho"]
}
mat_theta <- tmp
for(i in 1:27){
  mat_theta[tmp == i] <- canonicalPoints[i,"theta"]
}

# resolve NA's
rho_avg_list <- lapply(1:max(clone_assign), function(k){
  cell_idx <- which(clone_assign == k)
  vec <- matrixStats::colMeans2(mat_rho[cell_idx,], na.rm = T)
  vec[is.na(vec)] <- 1
  vec
})
theta_avg_list <- lapply(1:max(clone_assign), function(k){
  cell_idx <- which(clone_assign == k)
  vec <- matrixStats::colMeans2(mat_theta[cell_idx,],na.rm = T)
  vec[is.na(vec)] <- 0.5
  vec
})

for(k in 1:max(clone_assign)){
  cell_idx <- which(clone_assign == k)
  for(j in 1:ncol(mat_rho)){
    mat_rho[cell_idx[which(is.na(mat_rho[cell_idx,j]))], j] <- rho_avg_list[[k]][j]
  }
}

for(k in 1:max(clone_assign)){
  cell_idx <- which(clone_assign == k)
  for(j in 1:ncol(mat_theta)){
    mat_theta[cell_idx[which(is.na(mat_theta[cell_idx,j]))], j] <- theta_avg_list[[k]][j]
  }
}

# rescale the matrix
mat_rho <- mat_rho - 1
mat_rho <- mat_rho/svd(mat_rho)$d[1]*100
mat_theta <- mat_rho - 0.5
mat_theta <- mat_theta/svd(mat_theta)$d[1]*100
mat_2 <- cbind(mat_rho, mat_theta)
colnames(mat_2)[1:ncol(mat_rho)] <- paste0("rho-",colnames(mat_2)[1:ncol(mat_rho)])
colnames(mat_2)[ncol(mat_rho)+(1:ncol(mat_theta))] <- paste0("theta-",colnames(mat_2)[ncol(mat_rho)+(1:ncol(mat_theta))])

sd_vec <- matrixStats::colSds(mat_2)
if(any(sd_vec < 1e-6)){
  idx <- which(sd_vec < 1e-6)
  mat_2 <- mat_2[,-idx]
}

######################

all(colnames(mat) == names(clone_assign))
SNU <- Seurat::CreateSeuratObject(counts = mat,
                                  assay = "ATAC")
SNU[["clone"]] <- clone_assign
set.seed(10)
SNU <- Signac::RunTFIDF(SNU)
SNU <- Signac::FindTopFeatures(SNU, min.cutoff = 'q5')
SNU <- Signac::RunSVD(SNU, features = Seurat::VariableFeatures(object = SNU))
set.seed(10)
SNU <- Seurat::RunUMAP(SNU, reduction = 'lsi', dims = 1:10, 
                       reduction.name = "umap.atac", 
                       reduction.key = "atacUMAP_")

SNU[["CNA"]] <- Seurat::CreateAssayObject(counts = t(mat_2))
set.seed(10)
svd_res <- irlba::irlba(mat_2, nv = 50)
dimred <- multiomicCCA:::.mult_mat_vec(svd_res$u, svd_res$d)
rownames(dimred) <- colnames(SNU)
colnames(dimred) <- paste0("pca_", 1:ncol(dimred))
SNU[["pca"]] <- Seurat::CreateDimReducObject(embedding = dimred, 
                                             key = "pca_", assay = "CNA")
set.seed(10)
umap_res <- Seurat::RunUMAP(dimred)
umap_res <- umap_res@cell.embeddings
rownames(umap_res) <- colnames(SNU)
SNU[["umap.cna"]] <- Seurat::CreateDimReducObject(embedding = umap_res, 
                                                  key = "cnaUMAP_", assay = "CNA")

set.seed(10)
SNU <- Seurat::FindMultiModalNeighbors(SNU, 
                                       reduction.list = list("pca", "lsi"), 
                                       dims.list = list(1:50, 2:50))
set.seed(10)
SNU <- Seurat::RunUMAP(SNU, 
                       nn.name = "weighted.nn", 
                       reduction.name = "wnn.umap", 
                       reduction.key = "wnnUMAP_")

save(SNU, file = "../../../../out/Writeup14e/Writeup14e_SNU_preprocessed.RData")



