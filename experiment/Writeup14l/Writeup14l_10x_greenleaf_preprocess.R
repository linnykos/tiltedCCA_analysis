rm(list=ls())
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)


set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()


obj1 <- readRDS("~/nzhanglab/data/GSE162170_cortical_multiome/GSE162170_multiome_rna_counts.rds")
obj1 <- Matrix::Matrix(as.matrix(obj1), sparse = T)
obj2 <- readRDS("~/nzhanglab/data/GSE162170_cortical_multiome/GSE162170_multiome_atac_counts.rds")
obj2 <- Matrix::Matrix(as.matrix(obj2), sparse = T)

edb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
gene_name <- AnnotationDbi::mapIds(edb,
                                   keys = rownames(obj1),
                                   keytype = "GENEID", column = "GENENAME", multiVals = "first")
names(gene_name) <- NULL
na_idx <- which(is.na(gene_name))
if(length(na_idx) > 0){gene_name <- gene_name[-na_idx]; obj1 <- obj1[-na_idx,,drop=F]}
rownames(obj1) <- gene_name

colnames(obj2) <- colnames(obj1)
rownames(obj2) <- paste0("peak-", 1:nrow(obj2))

greenleaf <- Seurat::CreateSeuratObject(counts = obj1)
greenleaf[["ATAC"]] <- Seurat::CreateAssayObject(counts = obj2)

set.seed(10)
Seurat::DefaultAssay(greenleaf) <- "RNA"
greenleaf <- Seurat::SCTransform(greenleaf)
greenleaf <- Seurat::FindVariableFeatures(greenleaf)
greenleaf <- Seurat::RunPCA(greenleaf, verbose = FALSE)
set.seed(10)
greenleaf <- Seurat::RunUMAP(greenleaf, dims = 1:50)

####################

set.seed(10)
DefaultAssay(greenleaf) <- "ATAC"
greenleaf <- Signac::RunTFIDF(greenleaf)
greenleaf <-  Signac::FindTopFeatures(greenleaf, min.cutoff="q10")
greenleaf <-  Signac::RunSVD(greenleaf)  # question: is this svd run with only the top variable features?
set.seed(10)
greenleaf <- Seurat::RunUMAP(greenleaf, reduction="lsi", dims=2:50,reduction.name="umap.atac", 
                             reduction.key="atacUMAP_")

#####################

set.seed(10)
greenleaf <- Seurat::FindMultiModalNeighbors(greenleaf, reduction.list = list("pca", "lsi"), 
                                             dims.list = list(1:50, 2:50))
set.seed(10)
greenleaf <- Seurat::RunUMAP(greenleaf, nn.name = "weighted.nn", reduction.name = "umap.wnn", 
                             reduction.key = "wnnUMAP_")


cell_metadata <- readr::read_delim("~/nzhanglab/data/GSE162170_cortical_multiome/GSE162170_multiome_cell_metadata.txt.gz")
cluster_metadata <- readr::read_delim("~/nzhanglab/data/GSE162170_cortical_multiome/GSE162170_multiome_cluster_names.txt.gz")
cluster_metadata <- as.data.frame(cluster_metadata)
cluster_metadata <- cluster_metadata[which(cluster_metadata$Assay == "Multiome RNA"),]

all(cell_metadata$Cell.ID == colnames(greenleaf))
celltype_vec <- as.character(cell_metadata$seurat_clusters)
for(i in 1:nrow(cluster_metadata)){
  idx <- which(celltype_vec == cluster_metadata[i,"Cluster.ID"])
  celltype_vec[idx] <- cluster_metadata[i,"Cluster.Name"]
}
greenleaf$celltype <- celltype_vec
greenleaf$Sample.Age <- dplyr::pull(cell_metadata, "Sample.Age")
greenleaf$Sample.Batch <- dplyr::pull(cell_metadata, "Sample.Batch")

Seurat::DefaultAssay(greenleaf) <- "SCT"
greenleaf[["percent.mt"]] <- Seurat::PercentageFeatureSet(object = greenleaf, pattern = "^MT-")
greenleaf[["percent.rb"]] <- Seurat::PercentageFeatureSet(object = greenleaf, pattern = "^RPS")
greenleaf <- Seurat::CellCycleScoring(greenleaf, 
                                      g2m.features = cc.genes$g2m.genes, 
                                      s.features = cc.genes$s.genes)

#####################

save(greenleaf, date_of_run, session_info,
     file = "../../../../out/Writeup14l/Writeup14l_10x_greenleaf_preprocess.RData")

###############

plot1 <- Seurat::VlnPlot(greenleaf, features = c("nCount_ATAC", "nCount_RNA",
                                                 "percent.mt", "percent.rb"), ncol = 4,
                         pt.size = 0) + Seurat::NoLegend()
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14l/Writeup14l_10x_greenleaf_quality.png"),
                plot1, device = "png", width = 6.5, height = 5, units = "in")

png("../../../../out/figures/Writeup14l/Writeup14l_10x_greenleaf_seurat_umap.png", 
    height = 1500, width = 4500, units = "px", res = 300)
p1 <- Seurat::DimPlot(greenleaf, reduction = "umap", group.by = "celltype", 
                      label = TRUE, label.size = 2.5, repel = TRUE) + ggplot2::ggtitle("Greenleaf: RNA")
p2 <- Seurat::DimPlot(greenleaf, reduction = "umap.atac", group.by = "celltype", 
                      label = TRUE, label.size = 2.5, repel = TRUE) + ggplot2::ggtitle("Greenleaf: ATAC")
p3 <- Seurat::DimPlot(greenleaf, reduction = "umap.wnn", group.by = "celltype", 
                      label = TRUE, label.size = 2.5, repel = TRUE) + ggplot2::ggtitle("Greenleaf: WNN")
p1 + p2 + p3 & Seurat::NoLegend() & ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
graphics.off()

############################################################

library(dplyr)

# look at the gene activity matrix
Sys.setenv("VROOM_CONNECTION_SIZE" = 200000000)
mat <- readr::read_delim("~/nzhanglab/data/GSE162170_cortical_multiome/GSE162170_multiome_atac_gene_activities.tsv.gz",
                         delim = "\t")
## fix egregious error: the first column stores the gene names,
## but this makes all the columns offset by one
colname_vec <- colnames(mat)
rowname_vec <- dplyr::pull(mat, colname_vec[1])
mat <- mat[,-1]
last_vec <- dplyr::pull(mat, colname_vec[length(colname_vec)])
mat <- mat[,-ncol(mat)]
names(mat) <- colname_vec[1:ncol(mat)]
last_mat <- t(sapply(last_vec, function(x){
  as.numeric(strsplit(x, split = "\t")[[1]])
}))
colnames(last_mat) <- colname_vec[ncol(mat)+(1:2)]
rownames(last_mat) <- NULL
mat2 <- tibble::add_column(mat, last_mat[,1])
mat2 <- tibble::add_column(mat2, last_mat[,2])
names(mat2) <- colname_vec

colname_vec <- colnames(mat2)
sparse_mat_list <- lapply(1:length(colname_vec), function(i){
  if(i %% 100 == 0) print(paste0(i, " out of ", length(colname_vec)))
  Matrix::Matrix(dplyr::pull(mat2, colname_vec[i]), sparse = T)
})

len_vec <- sapply(sparse_mat_list, function(x){length(x@x)})
total_len <- sum(len_vec)
i_vec <- unlist(lapply(sparse_mat_list, function(x){x@i}))
j_vec <- unlist(lapply(1:length(sparse_mat_list), function(j){
  rep(j, length(sparse_mat_list[[j]]@i))
}))
x_vec <- unlist(lapply(sparse_mat_list, function(x){x@x}))

mat3 <- Matrix::sparseMatrix(i = i_vec+1,
                             j = j_vec,
                             x = x_vec,
                             dims = c(length(rowname_vec), length(colname_vec)))
rownames(mat3) <- rowname_vec
colnames(mat3) <- colname_vec

mat <- mat3
save(mat,
     file = "~/nzhanglab/data/GSE162170_cortical_multiome/GSE162170_multiome_atac_gene_activities.RData")
