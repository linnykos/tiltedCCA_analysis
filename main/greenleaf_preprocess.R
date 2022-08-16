rm(list=ls())

library(Seurat); library(Signac)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

create_seurat_object <- function(file_prefix, file_folder, file_suffix){
  tmp <- Seurat::Read10X_h5(paste0(file_prefix, file_folder, file_suffix))
  
  seurat_obj <- Seurat::CreateSeuratObject(
    counts = tmp[["Gene Expression"]],
    assay = "RNA"
  )
  
  seurat_obj
}

##########################

file_prefix <- "~/nzhanglab/data/GSE162170_cortical_multiome/raw-data/"
file_suffix <- "/filtered_feature_bc_matrix.h5"
file_folders <- c("hft_ctx_w21_dc1r3_r1", "hft_ctx_w21_dc2r2_r1", 
                  "hft_ctx_w21_dc2r2_r2")
name_vec <- c("hft_ctx_w21_dc1r3_r1", "hft_ctx_w21_dc2r2_r1", "hft_ctx_w21_dc2r2_r2")

annotation <- Signac::GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
ensembldb::seqlevelsStyle(annotation) <- "UCSC"
Signac::genome(annotation) <- "hg38"

seurat_list <- lapply(1:length(file_folders), function(i){
  print(i)
  file_folder <- file_folders[i]
  seurat_obj <- create_seurat_object(file_prefix = file_prefix,
                                     file_folder = file_folder,
                                     file_suffix = file_suffix)
  seurat_obj$dataset <- name_vec[i]
  print(dim(seurat_obj))
  seurat_obj
})

greenleaf <- merge(seurat_list[[1]], y = c(seurat_list[[2]], seurat_list[[3]]), 
                   add.cell.ids = name_vec, 
                   project = "hft", merge.data = T)

###################

# read in peak sets
peaks_list <- lapply(file_folders, function(file_folder){
  read.table(
    file = paste0(file_prefix, file_folder, "/atac_peaks.bed"),
    col.names = c("chr", "start", "end")
  )
})

# convert to genomic ranges
gr_list <- lapply(peaks_list, function(peaks_data){
  GenomicRanges::makeGRangesFromDataFrame(peaks_data)
})

# Create a unified set of peaks to quantify in each dataset
combined_peaks <- reduce(x = c(gr_list[[1]], gr_list[[2]], gr_list[[3]]))

# Filter out bad peaks based on length
peak_widths <- IRanges::width(combined_peaks)
combined_peaks <- combined_peaks[peak_widths < 10000 & peak_widths > 20]

# create fragment objects
frags_list <- lapply(1:length(file_folders), function(i){
  file_folder <- file_folders[i]
  Signac::CreateFragmentObject(
    path = paste0(file_prefix, file_folder, "/atac_fragments.tsv.gz"),
    cells = colnames(seurat_list[[i]])
  )
})

atac_count_list <- lapply(1:length(file_folders), function(i){
  Signac::FeatureMatrix(
    fragments = frags_list[[i]],
    features = combined_peaks,
    cells = colnames(seurat_list[[i]])
  )
})

seurat_atac_list <- lapply(1:length(file_folders), function(i){
  assay_obj <- Signac::CreateChromatinAssay(counts = atac_count_list[[i]],
                                            fragments = frags_list[[i]],
                                            sep = c("-", "-"),
                                            annotation = annotation)
  seurat_obj <- Seurat::CreateSeuratObject(assay_obj, assay = "ATAC")
  seurat_obj$dataset <- name_vec[i]
  seurat_obj
})

greenleaf_atac <- merge(seurat_atac_list[[1]], y = c(seurat_atac_list[[2]], seurat_atac_list[[3]]), 
                        add.cell.ids = name_vec, 
                        project = "hft", merge.data = T)

greenleaf[["ATAC"]] <- greenleaf_atac[["ATAC"]]
##########################

cell_metadata <- readr::read_delim("~/nzhanglab/data/GSE162170_cortical_multiome/GSE162170_multiome_cell_metadata.txt.gz")
cluster_metadata <- readr::read_delim("~/nzhanglab/data/GSE162170_cortical_multiome/GSE162170_multiome_cluster_names.txt.gz")
cluster_metadata <- as.data.frame(cluster_metadata)
cluster_metadata <- cluster_metadata[which(cluster_metadata$Assay == "Multiome RNA"),]
cell_metadata$Cell.ID <- paste0(cell_metadata$Cell.ID, "-1")
all(cell_metadata$Cell.ID %in% colnames(greenleaf))
keep_vec <- rep(0, ncol(greenleaf))
keep_vec[which(colnames(greenleaf) %in% cell_metadata$Cell.ID)] <- 1
greenleaf$keep <- keep_vec
greenleaf <- subset(greenleaf, keep == 1)
all(cell_metadata$Cell.ID == colnames(greenleaf))

celltype_vec <- as.character(cell_metadata$seurat_clusters)
for(i in 1:nrow(cluster_metadata)){
  idx <- which(celltype_vec == cluster_metadata[i,"Cluster.ID"])
  celltype_vec[idx] <- cluster_metadata[i,"Cluster.Name"]
}

greenleaf$celltype <- celltype_vec
greenleaf$Sample.Age <- dplyr::pull(cell_metadata, "Sample.Age")
greenleaf$Sample.Batch <- dplyr::pull(cell_metadata, "Sample.Batch")

###################

set.seed(10)
Seurat::DefaultAssay(greenleaf) <- "RNA"
greenleaf <- Seurat::SCTransform(greenleaf)
greenleaf <- Seurat::FindVariableFeatures(greenleaf)
greenleaf <- Seurat::RunPCA(greenleaf, verbose = FALSE)
set.seed(10)
greenleaf <- Seurat::RunUMAP(greenleaf, dims = 1:50)

set.seed(10)
DefaultAssay(greenleaf) <- "ATAC"
greenleaf <- Signac::RunTFIDF(greenleaf)
greenleaf <-  Signac::FindTopFeatures(greenleaf, min.cutoff="q10")
greenleaf <-  Signac::RunSVD(greenleaf)  
set.seed(10)
greenleaf <- Seurat::RunUMAP(greenleaf, reduction="lsi", 
                             dims=2:50, reduction.name="umap.atac", 
                             reduction.key="atacUMAP_")

####################

set.seed(10)
greenleaf <- Seurat::FindMultiModalNeighbors(greenleaf, reduction.list = list("pca", "lsi"), 
                                             dims.list = list(1:50, 2:50))
set.seed(10)
greenleaf <- Seurat::RunUMAP(greenleaf, nn.name = "weighted.nn", reduction.name = "umap.wnn", 
                             reduction.key = "wnnUMAP_")

Seurat::DefaultAssay(greenleaf) <- "SCT"
greenleaf[["percent.mt"]] <- Seurat::PercentageFeatureSet(object = greenleaf, pattern = "^MT-")
greenleaf[["percent.rb"]] <- Seurat::PercentageFeatureSet(object = greenleaf, pattern = "^RPS")
greenleaf <- Seurat::CellCycleScoring(greenleaf, 
                                      g2m.features = cc.genes$g2m.genes, 
                                      s.features = cc.genes$s.genes)

####################

# load in gene activity
load("~/nzhanglab/data/GSE162170_cortical_multiome/GSE162170_multiome_atac_gene_activities.RData")

name_vec <- colnames(greenleaf)
name_vec <- sapply(name_vec, function(x){strsplit(x, split = "-")[[1]][1]})
names(name_vec) <- NULL

mat <- mat[,name_vec]
colnames(mat) <- colnames(greenleaf)
rownames(mat) <- paste0("ATAC-", rownames(mat))

greenleaf[["geneActivity"]] <- Seurat::CreateAssayObject(counts = mat)
Seurat::DefaultAssay(greenleaf) <- "geneActivity"
greenleaf <- Seurat::NormalizeData(greenleaf)
greenleaf[["geneActivity"]]@var.features <- paste0("ATAC-", greenleaf[["SCT"]]@var.features)
greenleaf <- Seurat::ScaleData(greenleaf)
greenleaf <- Seurat::RunPCA(greenleaf, verbose = FALSE,
                            reduction.name = "pcaGeneActivity",
                            reduction.key = "PCAGeneActivity_")
set.seed(10)
greenleaf <- Seurat::RunUMAP(greenleaf, reduction="pcaGeneActivity", 
                             dims=2:50, reduction.name="umap.geneActivity", 
                             reduction.key = "geneActivityUMAP_")

###################

plot1 <- Seurat::DimPlot(greenleaf, reduction = "umap",
                         group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nRNA UMAP"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_greenleaf_rna-umap.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot2 <- Seurat::DimPlot(greenleaf, reduction = "umap.atac",
                         group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot2 <- plot2 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nATAC UMAP"))
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_greenleaf_atac-umap.png"),
                plot2, device = "png", width = 6, height = 5, units = "in")

plot3 <- Seurat::DimPlot(greenleaf, reduction = "umap.wnn",
                         group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot3 <- plot3 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nWNN UMAP"))
plot3 <- plot3 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_greenleaf_wnn-umap.png"),
                plot3, device = "png", width = 6, height = 5, units = "in")

plot4 <- Seurat::DimPlot(greenleaf, reduction = "umap.geneActivity",
                         group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot4 <- plot4 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nGene Activity UMAP"))
plot4 <- plot4 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_greenleaf_geneActivity-umap.png"),
                plot4, device = "png", width = 6, height = 5, units = "in")

###################

date_of_run <- Sys.time()
session_info <- devtools::session_info()

save(greenleaf, date_of_run, session_info,
     file = "../../../out/main/10x_greenleaf_preprocessed.RData")

#####################

Seurat::DefaultAssay(greenleaf) <- "ATAC"
set.seed(10)
greenleaf <- Seurat::FindNeighbors(greenleaf, dims = 1:50, reduction = "lsi")
set.seed(10)
greenleaf <- Seurat::FindClusters(greenleaf, resolution = 2)

plot1 <- Seurat::DimPlot(greenleaf, reduction = "umap.atac",
                         group.by = "seurat_clusters", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nATAC clusters"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_greenleaf_atac-clusters.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

for(i in 13:15){
  plot1 <- Seurat::DimPlot(greenleaf, reduction = "umap.atac",
                           cells.highlight = names(greenleaf$seurat_clusters[greenleaf$seurat_clusters == i]))
  plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nATAC clusters (Cluster ", i, ")"))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_greenleaf_atac-clusters_cluster", i, ".png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
}

######################################

# remove unwanted celltypes
greenleaf_tmp <- greenleaf
keep_vec <- rep(1, ncol(greenleaf_tmp))
keep_vec[greenleaf_tmp$celltype %in% c("EC/Peric.", "IN1", "IN2", "IN3", "mGPC/OPC", "SP")] <- 0
greenleaf_tmp$keep <- keep_vec
greenleaf_tmp <- subset(greenleaf_tmp, keep == 1)
table(greenleaf_tmp$seurat_clusters)

# remove unwanted clusters
keep_vec <- rep(1, ncol(greenleaf_tmp))
keep_vec[greenleaf_tmp$seurat_clusters %in% c("29", "22", "25", "9", "20", "11", "3", "6", "16", "27", "26", "30", "15")] <- 0
greenleaf_tmp$keep <- keep_vec
greenleaf_tmp <- subset(greenleaf_tmp, keep == 1)
table(greenleaf_tmp$seurat_clusters)

# ad-hoc merge some clusters for Slingshot to work better
greenleaf_tmp$seurat_clusters[greenleaf_tmp$seurat_clusters == "17"] <- "0"
greenleaf_tmp$seurat_clusters[greenleaf_tmp$seurat_clusters == "12"] <- "2"
greenleaf_tmp$seurat_clusters[greenleaf_tmp$seurat_clusters == "19"] <- "7"
greenleaf_tmp$seurat_clusters[greenleaf_tmp$seurat_clusters == "24"] <- "13"
greenleaf_tmp$seurat_clusters[greenleaf_tmp$seurat_clusters == "28"] <- "21"
greenleaf_tmp$seurat_clusters[greenleaf_tmp$seurat_clusters == "8"] <- "5"
table(greenleaf_tmp$seurat_clusters)

set.seed(10)
greenleaf_tmp <- Signac::RunSVD(greenleaf_tmp)  
set.seed(10)
greenleaf_tmp <- Seurat::FindNeighbors(greenleaf_tmp, dims = 2:50, reduction = "lsi")

dimmat_x <- greenleaf_tmp[["lsi"]]@cell.embeddings[,2:50]
snn_graph <- greenleaf_tmp@graphs$ATAC_nn

.nonzero_col <- function(mat, col_idx, bool_value){
  stopifnot(inherits(mat, "dgCMatrix"), col_idx %% 1 == 0,
            col_idx > 0, col_idx <= ncol(mat))
  
  val1 <- mat@p[col_idx]
  val2 <- mat@p[col_idx+1]
  
  if(val1 == val2) return(numeric(0))
  if(bool_value){
    # return the value
    mat@x[(val1+1):val2]
  } else {
    # return the row index
    mat@i[(val1+1):val2]+1
  }
}

# first get the MST, and then the curves
# If the lineages for the MST is undesireable, hard-set them yourself
# see 
# see https://github.com/LTLA/TrajectoryUtils/blob/master/R/createClusterMST.R (mainly calls to .create_mnn_distance_matrix and .estimate_edge_confidence)
set.seed(10)
slingshot_lin <- slingshot::getLineages(data = dimmat_x, 
                                        clusterLabels = as.character(greenleaf_tmp$seurat_clusters), 
                                        start.clus = c("10"), 
                                        end.clus = c("23", "13", "21"))
slingshot_lin@metadata$lineages
slingshot_curve <- slingshot::getCurves(slingshot_lin)

# extract pseudotime
pseudotime_mat <- slingshot_curve@assays@data$pseudotime
set.seed(10)
for(i in 1:3){
  # remove pseudotimes for cells not in the correct lineage
  idx <- greenleaf_tmp$seurat_clusters
  pseudotime_mat[which(!greenleaf_tmp$seurat_clusters %in% slingshot_lin@metadata$lineages[[i]]),i] <- NA
  idx <- which(!is.na(pseudotime_mat[,i]))
  pseudotime_mat[idx,i] <- rank(pseudotime_mat[idx,i])
  pseudotime_mat_tmp <- pseudotime_mat
  for(k in idx){
    neighbors <- .nonzero_col(snn_graph, col_idx = k, bool_value = F)
    neighbors <- unique(c(intersect(neighbors, idx), k))
    pseudotime_mat_tmp[k,i] <- mean(pseudotime_mat[neighbors,i], na.rm = T)
  }
  pseudotime_mat <- pseudotime_mat_tmp
  pseudotime_mat[idx,i] <- rank(pseudotime_mat[idx,i], ties.method = "random")
}
# compute a weighted rank
lineage_length <- sapply(slingshot_lin@metadata$lineages, function(lin){
  length(which(greenleaf_tmp$seurat_clusters %in% lin))
})
pseudotime_vec <- apply(pseudotime_mat, 1, function(x){
  vec <- sapply(1:3, function(i){x[i]/lineage_length[i]})
  mean(vec, na.rm = T)
})
pseudotime_vec_full <- rep(NA, ncol(greenleaf))
names(pseudotime_vec_full) <- colnames(greenleaf)
pseudotime_vec_full[names(pseudotime_vec)] <- pseudotime_vec
greenleaf$pseudotime <- pseudotime_vec_full

# add metadata
for(i in 1:3){
  traj_vec <- rep(0, ncol(greenleaf))
  names(traj_vec) <- colnames(greenleaf)
  idx <- which(!is.na(pseudotime_mat[,i]))
  traj_vec[names(pseudotime_mat[idx,i])] <- 1
  greenleaf <- Seurat::AddMetaData(greenleaf, metadata = traj_vec, col.name = paste0("Lineage", i))
}

###############

# consensus pca
Seurat::DefaultAssay(greenleaf) <- "SCT"
mat_1 <- Matrix::t(greenleaf[["SCT"]]@data[Seurat::VariableFeatures(object = greenleaf),])
Seurat::DefaultAssay(greenleaf) <- "ATAC"
mat_2 <- Matrix::t(greenleaf[["ATAC"]]@data[Seurat::VariableFeatures(object = greenleaf),])

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

consensus_pca <- tiltedCCA:::consensus_pca(mat_1 = mat_1b, mat_2 = mat_2b,
                                           dims_1 = 1:50, dims_2 = 2:50,
                                           dims_consensus = 1:49,
                                           center_1 = T, center_2 = F,
                                           recenter_1 = F, recenter_2 = T,
                                           rescale_1 = F, rescale_2 = T,
                                           scale_1 = T, scale_2 = F,
                                           verbose = 1)
consensus_dimred <- consensus_pca$dimred_consensus
colnames(consensus_dimred) <- paste0("consensusPCA_", 1:ncol(consensus_dimred))
greenleaf[["consensusPCA"]] <- Seurat::CreateDimReducObject(consensus_dimred, 
                                                            assay = "SCT")

set.seed(10)
umap_res <- Seurat::RunUMAP(consensus_pca$dimred_consensus)
umap_mat <- umap_res@cell.embeddings
rownames(umap_mat) <- colnames(greenleaf)
colnames(umap_mat) <- paste0("consensusUMAP_", 1:ncol(umap_mat))
greenleaf[["consensusUMAP"]] <- Seurat::CreateDimReducObject(umap_mat, 
                                                             assay = "SCT")

##########

save(greenleaf, date_of_run, session_info,
     file = "../../../out/main/10x_greenleaf_preprocessed.RData")

plot1 <- Seurat::FeaturePlot(greenleaf, reduction = "umap.atac",
                             features = "pseudotime")
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC): ATAC UMAP\nSlingshot's pseudotime via ATAC"))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_greenleaf_atac-pseudotime.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <- Seurat::FeaturePlot(greenleaf, reduction = "umap",
                             features = "pseudotime")
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC): RNA UMAP\nSlingshot's pseudotime via ATAC"))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_greenleaf_rna-pseudotime.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot3 <- Seurat::FeaturePlot(greenleaf, reduction = "umap.wnn",
                             features = "pseudotime")
plot3 <- plot3 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC): WNN UMAP\nSlingshot's pseudotime via ATAC"))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_greenleaf_wnn-pseudotime.png"),
                plot3, device = "png", width = 6, height = 5, units = "in")

plot4 <- Seurat::FeaturePlot(greenleaf, reduction = "umap.geneActivity",
                             features = "pseudotime")
plot4 <- plot4 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC): Gene Activity UMAP\nSlingshot's pseudotime via ATAC"))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_greenleaf_geneActivity-pseudotime.png"),
                plot4, device = "png", width = 6, height = 5, units = "in")

png("../../../out/figures/main/10x_greenleaf_pseudotime_histogram.png",
    height = 1500, width = 2500, res = 300, units = "px")
hist(greenleaf$pseudotime, main = "Slingshot pseudotime", xlab = "Pseudotime", breaks = 50)
graphics.off()


plot1 <- Seurat::DimPlot(greenleaf, reduction = "consensusUMAP",
                         group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nConsensus PCA's UMAP"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_greenleaf_consensusPCA-umap.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")


######################

greenleaf$logumi <- log1p(greenleaf$nCount_RNA)
plot1 <- Seurat::VlnPlot(greenleaf, features = "logumi", 
                         group.by = 'celltype', 
                         sort = "decreasing") + Seurat::NoLegend()
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_greenleaf_ncount_violin-by-celltype.png"),
                plot1, device = "png", width = 10, height = 3, units = "in",
                dpi = 500)
