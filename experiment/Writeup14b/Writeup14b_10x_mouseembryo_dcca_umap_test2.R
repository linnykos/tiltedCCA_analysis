rm(list=ls()); set.seed(10)

library(Seurat); library(Signac); library(EnsDb.Hsapiens.v86); 
library(BSgenome.Hsapiens.UCSC.hg38); library(GenomeInfoDb)
library(dplyr); library(ggplot2); library(multiomicCCA)

load("../../../../out/Writeup14/Writeup14_10x_mouseembryo_preprocess.RData")
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

##########

set.seed(10)
obj <- multiomicCCA:::.prepare_embeddings_singleton(dcca_res$common_score, dcca_res$distinct_score_1,
                                             dcca_res$svd_1, add_noise = F)
set.seed(10)
seurat_umap <- Seurat::RunUMAP(obj$common, metric = "euclidean")

mbrain2 <- Seurat::CreateSeuratObject(counts = t(mat_1))
mbrain2[["label_Savercat"]] <- metadata$label_Savercat[cell_idx]
mbrain2[["asdf"]] <- Seurat::CreateDimReducObject(embedding = seurat_umap@cell.embeddings, 
                                                  key = "UMAP", assay = "RNA")
plot1 <- Seurat::DimPlot(mbrain2, reduction = "asdf", label = TRUE,
                         group.by = "label_Savercat",
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_umap_original_test.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

##################

# now try recreating it via the graph
set.seed(10)
rann_res <- RANN::nn2(obj$common, k = 30)
rann_res$nn.idx <- rann_res$nn.idx[,-1]
rann_res$nn.dist <- rann_res$nn.dist[,-1]

rann_res$nn.idx <- lapply(1:nrow(rann_res$nn.idx), function(i){
  rann_res$nn.idx[i,]
})
rann_res$nn.dist <- lapply(1:nrow(rann_res$nn.dist), function(i){
  rann_res$nn.dist[i,]
})

j_vec <- unlist(rann_res$nn.idx)
i_vec <- unlist(lapply(1:length(rann_res$nn.idx), function(i){
  rep(i, length(rann_res$nn.idx[[i]]))
}))
x_vec <- unlist(rann_res$nn.dist)
dist_mat <- Matrix::sparseMatrix(i = i_vec, j = j_vec, x = x_vec)

rownames(dist_mat) <- rownames(mat_1)
colnames(dist_mat) <- rownames(mat_1)
g_obj <- SeuratObject::as.Graph(dist_mat)

# what if we use euclidean
set.seed(10)
seurat_umap2 <- Seurat::RunUMAP(g_obj, metric = "euclidean")
mbrain2[["asdf"]] <- Seurat::CreateDimReducObject(embedding = seurat_umap2@cell.embeddings, 
                                                  key = "UMAP", assay = "RNA")
plot1 <- Seurat::DimPlot(mbrain2, reduction = "asdf", label = TRUE,
                         group.by = "label_Savercat",
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_umap_original_test2.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

# what if we use cosine
set.seed(10)
seurat_umap2 <- Seurat::RunUMAP(g_obj, metric = "cosine")
mbrain2[["asdf"]] <- Seurat::CreateDimReducObject(embedding = seurat_umap2@cell.embeddings, 
                                                  key = "UMAP", assay = "RNA")
plot1 <- Seurat::DimPlot(mbrain2, reduction = "asdf", label = TRUE,
                         group.by = "label_Savercat",
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_umap_original_test3.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")


############
# don't remove the diagonal
set.seed(10)
rann_res <- RANN::nn2(obj$common, k = 30)
rann_res$nn.idx <- rann_res$nn.idx[,-1]
rann_res$nn.dist <- rann_res$nn.dist[,-1]
rann_res$nn.idx <- lapply(1:nrow(rann_res$nn.idx), function(i){
  rann_res$nn.idx[i,]
})
rann_res$nn.dist <- lapply(1:nrow(rann_res$nn.dist), function(i){
  rann_res$nn.dist[i,]
})
j_vec <- unlist(rann_res$nn.idx)
i_vec <- unlist(lapply(1:length(rann_res$nn.idx), function(i){
  rep(i, length(rann_res$nn.idx[[i]]))
}))
x_vec <- unlist(rann_res$nn.dist)
dist_mat <- Matrix::sparseMatrix(i = i_vec, j = j_vec, x = x_vec)
rownames(dist_mat) <- rownames(mat_1)
colnames(dist_mat) <- rownames(mat_1)
g_obj <- SeuratObject::as.Graph(dist_mat)
set.seed(42)
reticulate::py_set_seed(42)
seurat_umap2 <- Seurat:::RunUMAP.Graph(g_obj, metric = "euclidean", seed.use = 10L)
mbrain2[["asdf"]] <- Seurat::CreateDimReducObject(embedding = seurat_umap2@cell.embeddings, 
                                                  key = "UMAP", assay = "RNA")
plot1 <- Seurat::DimPlot(mbrain2, reduction = "asdf", label = TRUE,
                         group.by = "label_Savercat",
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_umap_original_test4.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")
set.seed(42)
reticulate::py_set_seed(42)
seurat_umap2 <- Seurat:::RunUMAP.Graph(g_obj, metric = "euclidean", seed.use = 10L)
mbrain2[["asdf"]] <- Seurat::CreateDimReducObject(embedding = seurat_umap2@cell.embeddings, 
                                                  key = "UMAP", assay = "RNA")
plot1 <- Seurat::DimPlot(mbrain2, reduction = "asdf", label = TRUE,
                         group.by = "label_Savercat",
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_umap_original_test4b.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

#######################

# what if we called reticulate ourselves?
library(reticulate)
reticulate::use_virtualenv("~/.virtualenvs")
object <- g_obj
np <- reticulate::import("numpy", delay_load = TRUE)
sp <- reticulate::import("scipy", delay_load = TRUE)
sklearn <- reticulate::import("sklearn", delay_load = TRUE)
umap <- reticulate::import("umap", delay_load = TRUE)
diag(x = object) <- 0
data <- object
object <- sp$sparse$coo_matrix(arg1 = object)
umap.method = 'umap-learn'
n.components = 2L
metric = 'correlation'
n.epochs = 0L
learning.rate = 1
min.dist = 0.3
spread = 1
repulsion.strength = 1
negative.sample.rate = 5L
a = NULL
b = NULL
uwot.sgd = FALSE
seed.use = 42L
metric.kwds = NULL
verbose = TRUE
reduction.key = 'UMAP_'
ab.params <- umap$umap_$find_ab_params(spread = spread, min_dist = min.dist)
a <- a %||% ab.params[[1]]
b <- b %||% ab.params[[2]]
n.epochs <- n.epochs %||% 0L
#random.state <- sklearn$utils$check_random_state(seed = as.integer(x = 42))
random.state <- np$random$RandomState(42L)
umap.args <- list(
  data = object,
  graph = object,
  n_components = n.components,
  initial_alpha = learning.rate,
  a = a,
  b = b,
  gamma = repulsion.strength,
  negative_sample_rate = negative.sample.rate,
  n_epochs = as.integer(x = n.epochs),
  random_state = random.state,
  init = "spectral",
  metric = metric,
  metric_kwds = metric.kwds,
  verbose = verbose
)
umap.args <- c(umap.args, list(
  densmap = FALSE,
  densmap_kwds = NULL,
  output_dens = FALSE
))
embeddings <- do.call(what = umap$umap_$simplicial_set_embedding, args = umap.args)
if (length(x = embeddings) == 2) {
  embeddings <- embeddings[[1]]
}
rownames(x = embeddings) <- colnames(x = data)
colnames(x = embeddings) <- paste0("UMAP_", 1:n.components)
# center the embeddings on zero
embeddings <- scale(x = embeddings, scale = FALSE)
umap <- CreateDimReducObject(
  embeddings = embeddings,
  key = reduction.key,
  assay = NULL,
  global = TRUE
)
mbrain2[["asdf"]] <- Seurat::CreateDimReducObject(embedding = umap@cell.embeddings, 
                                                  key = "UMAP", assay = "RNA")
plot1 <- Seurat::DimPlot(mbrain2, reduction = "asdf", label = TRUE,
                         group.by = "label_Savercat",
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_umap_original_test4.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

library(reticulate)
object <- g_obj
np <- reticulate::import("numpy", delay_load = TRUE)
sp <- reticulate::import("scipy", delay_load = TRUE)
sklearn <- reticulate::import("sklearn", delay_load = TRUE)
umap <- reticulate::import("umap", delay_load = TRUE)
diag(x = object) <- 0
data <- object
object <- sp$sparse$coo_matrix(arg1 = object)
umap.method = 'umap-learn'
n.components = 2L
metric = 'correlation'
n.epochs = 0L
learning.rate = 1
min.dist = 0.3
spread = 1
repulsion.strength = 1
negative.sample.rate = 5L
a = NULL
b = NULL
uwot.sgd = FALSE
seed.use = 42L
metric.kwds = NULL
verbose = TRUE
reduction.key = 'UMAP_'
ab.params <- umap$umap_$find_ab_params(spread = spread, min_dist = min.dist)
a <- a %||% ab.params[[1]]
b <- b %||% ab.params[[2]]
n.epochs <- n.epochs %||% 0L
#random.state <- sklearn$utils$check_random_state(seed = as.integer(x = 42))
random.state <- np$random$RandomState(42L)
umap.args <- list(
  data = object,
  graph = object,
  n_components = n.components,
  initial_alpha = learning.rate,
  a = a,
  b = b,
  gamma = repulsion.strength,
  negative_sample_rate = negative.sample.rate,
  n_epochs = as.integer(x = n.epochs),
  random_state = random.state,
  init = "spectral",
  metric = metric,
  metric_kwds = metric.kwds,
  verbose = verbose
)
umap.args <- c(umap.args, list(
  densmap = FALSE,
  densmap_kwds = NULL,
  output_dens = FALSE
))
embeddings <- do.call(what = umap$umap_$simplicial_set_embedding, args = umap.args)
if (length(x = embeddings) == 2) {
  embeddings <- embeddings[[1]]
}
rownames(x = embeddings) <- colnames(x = data)
colnames(x = embeddings) <- paste0("UMAP_", 1:n.components)
# center the embeddings on zero
embeddings <- scale(x = embeddings, scale = FALSE)
umap <- CreateDimReducObject(
  embeddings = embeddings,
  key = reduction.key,
  assay = NULL,
  global = TRUE
)
mbrain2[["asdf"]] <- Seurat::CreateDimReducObject(embedding = umap@cell.embeddings, 
                                                  key = "UMAP", assay = "RNA")
plot1 <- Seurat::DimPlot(mbrain2, reduction = "asdf", label = TRUE,
                         group.by = "label_Savercat",
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_umap_original_test4b.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

