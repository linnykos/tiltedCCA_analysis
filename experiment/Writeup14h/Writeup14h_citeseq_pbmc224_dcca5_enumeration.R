rm(list=ls())
load("../../../../out/Writeup14h/Writeup14h_citeseq_pbmc224_preprocessed.RData")
library(Seurat)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()
n <- ncol(pbmc)

membership_vec <- as.factor(pbmc@meta.data$celltype.l2)
n <- length(membership_vec)
sort(table(membership_vec))
set.seed(10)
idx <- multiomicCCA::construct_celltype_subsample(membership_vec, 
                                                  min_subsample_cell = 5000)
keep_vec <- rep(0, n)
names(keep_vec) <- rownames(pbmc@meta.data)
keep_vec[idx] <- 1
pbmc$keep <- keep_vec
pbmc <- subset(pbmc, keep == 1)
pbmc[["rna.umap"]] <- NULL
pbmc[["adt.umap"]] <- NULL
pbmc[["wnn.umap"]] <- NULL

set.seed(10)
pbmc <- Seurat::RunUMAP(pbmc, reduction = 'pca', dims = 1:40, assay = 'SCT',
                        reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
set.seed(10)
pbmc <- Seurat::RunUMAP(pbmc, reduction = 'apca', dims = 1:50, assay = 'ADT',
                        reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')

###############

mat_1 <- t(pbmc[["SCT"]]@scale.data)
mat_2 <- t(pbmc[["ADT"]]@scale.data)

mat_1b <- mat_1
sd_vec <- sparseMatrixStats::colSds(mat_1b)
if(any(sd_vec <= 1e-6)){
  mat_1b <- mat_1b[,-which(sd_vec <= 1e-6)]
}

mat_2b <- mat_2
sd_vec <- sparseMatrixStats::colSds(mat_2b)
if(any(sd_vec <= 1e-6)){
  mat_2b <- mat_2b[,-which(sd_vec <= 1e-6)]
}

#########################

rank_1 <- 40; rank_2 <- 50; nn <- 30
discretization_gridsize <- rev(seq(0, 1, length.out = 11))
dcca_list <- vector("list", length = length(discretization_gridsize))
for(i in 1:length(discretization_gridsize)){
  tilt <- discretization_gridsize[i]
  print(paste0("Starting tilt: ", tilt))
  
  set.seed(10)
  dcca_list[[i]] <- multiomicCCA::dcca_factor(mat_1b, mat_2b, 
                                              dims_1 = 1:rank_1, dims_2 = 1:rank_2,
                                              center_1 = F, center_2 = F,
                                              scale_1 = F, scale_2 = F,
                                              num_neigh = nn, 
                                              metacell_clustering_1 = NA,
                                              metacell_clustering_2 = NA,
                                              fix_tilt_perc = tilt, 
                                              enforce_boundary = F,
                                              verbose = T)
  
  save(dcca_list,
       file = "../../../../out/Writeup14h/Writeup14h_citeseq_pbmc224_dcca5_enumeration_tmp.RData")
}

save(pbmc, dcca_list, 
     rank_1, rank_2, nn, date_of_run, session_info,
     discretization_gridsize,
     file = "../../../../out/Writeup14h/Writeup14h_citeseq_pbmc224_dcca5_enumeration.RData")
