rm(list=ls())

library(Signac)
library(Seurat)
library(EnsDb.Mmusculus.v79)
library(GenomeInfoDb)
library(ggplot2)
library(patchwork)
library(biovizBase)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

# --------------------------------------------------
# Load in data and set up Seurat object
# --------------------------------------------------

# Create a Seurat object first with the RNA
counts <- Seurat::Read10X_h5("~/nzhanglab/data/10x_mouse_brain/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix.h5")
rownames(counts[["Gene Expression"]]) <- toupper(rownames(counts[["Gene Expression"]]))
metadata <- read.csv(
  file = "~/nzhanglab/data/10x_mouse_brain/e18_mouse_brain_fresh_5k_per_barcode_metrics.csv",
  header = TRUE,
  row.names = 1
)

mbrain <- Seurat::CreateSeuratObject(counts = counts[["Gene Expression"]],
                                     meta.data = metadata)


# Create an assay for the ATAC
# only use peaks in standard chromosomes.
grange.counts <- Signac::StringToGRanges(rownames(counts$Peaks), sep=c(":","-"))
grange.use <- GenomeInfoDb::seqnames(grange.counts) %in% GenomeInfoDb::standardChromosomes(grange.counts)
atac_counts <- counts[["Peaks"]][as.vector(grange.use),]

# needs the e18_mouse_brain_fresh_5k_atac_fragments.tsv.gz.tbi file also
chrom_assay <- Signac::CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = "~/nzhanglab/data/10x_mouse_brain/e18_mouse_brain_fresh_5k_atac_fragments.tsv.gz",
  min.cells = 10
)

# This step is VERY important! Set the annotations levels style so that
# "counts" and "annotations" both use the same name for chromosomes.
# that is, "chr1" and not "1".
annotations <- Signac::GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79) # takes a while
GenomeInfoDb::seqlevelsStyle(annotations) <- "UCSC" # chromosome 1 is called "chr1"
annotations
Signac::Annotation(chrom_assay) <- annotations

mbrain[["ATAC"]] <- chrom_assay # add to Seurat object.

# --------------------------------------------------
# Gene activities
# --------------------------------------------------

set.seed(10)
Seurat::DefaultAssay(mbrain) <- "ATAC"
gene_activities <- Signac::GeneActivity(mbrain,
                                        extend.upstream = 500,
                                        extend.downstream = 500)
rownames(gene_activities) <- toupper(rownames(gene_activities))
mbrain[["geneActivity"]] <- Seurat::CreateAssayObject(counts = gene_activities)

# --------------------------------------------------
# Get transferred cell labels from Jane.
# --------------------------------------------------

mbrain.jane <- readRDS("~/project/tiltedCCA/data/10x_mouseembryo/data_tenx_labels_Jane.rds")
mbrain$label_Savercat <- mbrain.jane$savercatLable

# --------------------------------------------------
# Preprocess RMA
# --------------------------------------------------

set.seed(10)
Seurat::DefaultAssay(mbrain) <- "RNA"
mbrain <- Seurat::SCTransform(mbrain)
mbrain <- Seurat::FindVariableFeatures(mbrain)

# --------------------------------------------------
# Preprocess Gene Activity
# --------------------------------------------------

set.seed(10)
Seurat::DefaultAssay(mbrain) <- "geneActivity"
mbrain <- Seurat::NormalizeData(mbrain)
var_features <- mbrain[["SCT"]]@var.features
mbrain[["geneActivity"]]@var.features <- intersect(var_features, rownames(mbrain))
mbrain <- Seurat::ScaleData(mbrain)

# --------------------------------------------------
# Preprocess ATAC
# --------------------------------------------------

set.seed(10)
DefaultAssay(mbrain) <- "ATAC"
mbrain <- Signac::RunTFIDF(mbrain)
mbrain <-  Signac::FindTopFeatures(mbrain, min.cutoff="q10")

# --------------------------------------------------
# Subselect only the relevant cells
# --------------------------------------------------

n <- ncol(mbrain)
keep_vec <- rep(0, n)
keep_vec[which(mbrain$label_Savercat %in% c("Oligodendrocyte", "Radial glia", 
                                            "Forebrain GABAergic", "Neuroblast", 
                                            "Glioblast", "Cortical or hippocampal glutamatergic"))] <- 1
mbrain$keep <- keep_vec
mbrain <- subset(mbrain, keep == 1)

# --------------------------------------------------
# RNA dimension-reduction
# --------------------------------------------------

set.seed(10)
DefaultAssay(mbrain) <- "SCT"
mbrain <- Seurat::RunPCA(mbrain, verbose = FALSE,
                         reduction.name = "pca",
                         reduction.key = "PCA_")
set.seed(10)
mbrain <- Seurat::RunUMAP(mbrain, reduction="pca", 
                          dims = 1:50, 
                          reduction.name="umap", 
                          reduction.key = "UMAP_")

# --------------------------------------------------
# Gene activity dimension-reduction
# --------------------------------------------------

set.seed(10)
DefaultAssay(mbrain) <- "geneActivity"
mbrain <- Seurat::RunPCA(mbrain, verbose = FALSE,
                         reduction.name = "pca.ga",
                         reduction.key = "geneActPCA_")
set.seed(10)
mbrain <- Seurat::RunUMAP(mbrain, reduction="pca.ga", 
                          dims = 2:50, 
                          reduction.name="umap.ga", 
                          reduction.key = "geneActUMAP_")

# --------------------------------------------------
# ATAC dimension-reduction
# --------------------------------------------------

set.seed(10)
DefaultAssay(mbrain) <- "ATAC"
mbrain <-  Signac::RunSVD(mbrain)  
set.seed(10)
mbrain <- Seurat::RunUMAP(mbrain, 
                          reduction="lsi", 
                          dims=2:50,
                          reduction.name="umap.atac", 
                          reduction.key="atacUMAP_")

# --------------------------------------------------
# Multi-modal analyses (WNN and Consensus PCA)
# --------------------------------------------------

set.seed(10)
mbrain <- Seurat::FindMultiModalNeighbors(mbrain, 
                                          reduction.list = list("pca", "lsi"), 
                                          dims.list = list(1:50, 2:50))
set.seed(10)
mbrain <- Seurat::RunUMAP(mbrain, 
                          nn.name = "weighted.nn", 
                          reduction.name = "umap.wnn", 
                          reduction.key = "wnnUMAP_")

Seurat::DefaultAssay(mbrain) <- "SCT"
mat_1 <- Matrix::t(mbrain[["SCT"]]@data[Seurat::VariableFeatures(object = mbrain),])
Seurat::DefaultAssay(mbrain) <- "ATAC"
mat_2 <- Matrix::t(mbrain[["ATAC"]]@data[Seurat::VariableFeatures(object = mbrain),])

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
mbrain[["consensusPCA"]] <- Seurat::CreateDimReducObject(consensus_dimred, 
                                                         assay = "SCT")

set.seed(10)
umap_res <- Seurat::RunUMAP(consensus_pca$dimred_consensus)
umap_mat <- umap_res@cell.embeddings
rownames(umap_mat) <- colnames(mbrain)
colnames(umap_mat) <- paste0("consensusUMAP_", 1:ncol(umap_mat))
mbrain[["consensusUMAP"]] <- Seurat::CreateDimReducObject(umap_mat, 
                                                          assay = "SCT")

# --------------------------------------------------
# Label slingshot pseudotime
# --------------------------------------------------

Seurat::DefaultAssay(mbrain) <- "ATAC"
set.seed(10)
mbrain <- Seurat::FindNeighbors(mbrain, dims = 2:50, reduction = "lsi")
set.seed(10)
mbrain <- Seurat::FindClusters(mbrain, resolution = 2)

# remove unwanted clusters
keep_vec <- rep(1, ncol(mbrain))
keep_vec[mbrain$seurat_clusters %in% c("5","15","19","6","18","13","14","16", "4")] <- 0
mbrain$keep <- keep_vec
mbrain_tmp <- subset(mbrain, keep == 1)
mbrain_tmp2 <- mbrain_tmp
table(mbrain_tmp$seurat_clusters)

# ad-hoc merge some clusters for Slingshot to work better
mbrain_tmp$seurat_clusters[mbrain_tmp$seurat_clusters == "1"] <- "8"
mbrain_tmp$seurat_clusters[mbrain_tmp$seurat_clusters == "2"] <- "8"
mbrain_tmp$seurat_clusters[mbrain_tmp$seurat_clusters == "20"] <- "9"
mbrain_tmp$seurat_clusters[mbrain_tmp$seurat_clusters == "11"] <- "0"

set.seed(10)
mbrain_tmp <-  Signac::RunSVD(mbrain_tmp)  
set.seed(10)
mbrain_tmp <- Seurat::FindNeighbors(mbrain_tmp, dims = 2:50, reduction = "lsi")

dimmat_x <- mbrain_tmp[["lsi"]]@cell.embeddings[,2:50]
snn_graph <- mbrain_tmp@graphs$ATAC_nn

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
                                        clusterLabels = as.character(mbrain_tmp$seurat_clusters), 
                                        start.clus = c("17"), 
                                        end.clus = c("21", "12", "8"))
slingshot_lin@metadata$lineages
slingshot_curve <- slingshot::getCurves(slingshot_lin)

# minor adjustment
.initial_curve_fit <- function(cluster_vec,
                               dimred, 
                               lineage_order){
  stopifnot(all(lineage_order %in% cluster_vec),
            length(cluster_vec) == nrow(dimred))
  t(sapply(lineage_order, function(cluster){
    idx <- which(cluster_vec == cluster)
    Matrix::colMeans(dimred[idx,,drop = F])
  }))
}

.extract_pseudotime <- function(dimred,
                                initial_fit,
                                stretch){ # default strecth=2
  pcurve <- princurve::project_to_curve(dimred,
                                        s = initial_fit,
                                        stretch = 2)
  pcurve$lambda
}

lineage1_initial_fit <- .initial_curve_fit(cluster_vec = as.character(mbrain_tmp2$seurat_clusters),
                                           dimred = dimmat_x,
                                           lineage_order = c("17", "9", "3", "0", "10", "8"))
lineage1_pseudotime <- .extract_pseudotime(dimred = dimmat_x[which(mbrain_tmp$seurat_clusters %in% slingshot_lin@metadata$lineages[[1]]),],
                                           initial_fit = lineage1_initial_fit,
                                           stretch = 2)
# sapply(unique(as.character(c("17", "9", "3", "11", "0", "10", "1", "2", "8"))), function(clust){
#   cellnames <- colnames(mbrain_tmp2)[which(mbrain_tmp2$seurat_clusters == clust)]
#   mean(lineage1_pseudotime[cellnames])
# })

# extract pseudotime
pseudotime_mat <- slingshot_curve@assays@data$pseudotime
num_lin <- length(slingshot_lin@metadata$lineages)
pseudotime_mat[names(lineage1_pseudotime),1] <- lineage1_pseudotime

for(i in 1:num_lin){
  # remove pseudotimes for cells not in the correct lineage
  pseudotime_mat[which(!mbrain_tmp$seurat_clusters %in% slingshot_lin@metadata$lineages[[i]]),i] <- NA
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
  length(which(mbrain_tmp$seurat_clusters %in% lin))
})
pseudotime_vec <- apply(pseudotime_mat, 1, function(x){
  vec <- sapply(1:num_lin, function(i){x[i]/lineage_length[i]})
  mean(vec, na.rm = T)
})

# store results
pseudotime_vec_full <- rep(NA, ncol(mbrain))
names(pseudotime_vec_full) <- colnames(mbrain)
pseudotime_vec_full[names(pseudotime_vec)] <- pseudotime_vec
mbrain$pseudotime <- pseudotime_vec_full
# table(table(pseudotime_vec))

# add metadata
# mbrain$Lineage1 <- NULL; mbrain$Lineage2 <- NULL; mbrain$Lineage3 <- NULL
lineage_celltypes <- list(c("Radial glia", "Neuroblast", "Cortical or hippocampal glutamatergic"),
                          c("Radial glia", "Neuroblast", "Forebrain GABAergic"),
                          c("Radial glia", "Glioblast", "Oligodendrocyte"))
for(i in 1:num_lin){
  traj_vec <- rep(0, ncol(mbrain))
  names(traj_vec) <- colnames(mbrain)
  
  cellnames1 <- colnames(mbrain_tmp)[which(mbrain_tmp$label_Savercat %in% lineage_celltypes[[i]])]
  cellnames2 <- colnames(mbrain_tmp)[which(mbrain_tmp$seurat_clusters %in% slingshot_lin@metadata$lineages[[i]])]
  
  cellnames <- intersect(cellnames1, cellnames2)
  traj_vec[cellnames] <- 1
  mbrain <- Seurat::AddMetaData(mbrain, metadata = traj_vec, col.name = paste0("Lineage", i))
}

########################################
########################################
########################################

save(mbrain, date_of_run, session_info,
     file = "../../../out/main/10x_mouseembryo_preprocessed.RData")

########################################
########################################
########################################

source("mouseembryo_colorPalette.R")

plot1 <- Seurat::DimPlot(mbrain, reduction = "umap",
                         group.by = "label_Savercat", label = TRUE,
                         cols = col_palette, 
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Mouse Embryo E18 (10x, RNA+ATAC)\nRNA UMAP"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_mouseembryo_rna-umap.png"),
                plot1, device = "png", width = 7, height = 5, units = "in")

plot2 <- Seurat::DimPlot(mbrain, reduction = "umap.ga",
                         group.by = "label_Savercat", label = TRUE,
                         cols = col_palette, 
                         repel = TRUE, label.size = 2.5)
plot2 <- plot2 + ggplot2::ggtitle(paste0("Mouse Embryo E18 (10x, RNA+ATAC)\nGene Activity UMAP"))
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_mouseembryo_geneActivity-umap.png"),
                plot2, device = "png", width = 7, height = 5, units = "in")

plot2 <- Seurat::DimPlot(mbrain, reduction = "umap.atac",
                         group.by = "label_Savercat", label = TRUE,
                         cols = col_palette, 
                         repel = TRUE, label.size = 2.5)
plot2 <- plot2 + ggplot2::ggtitle(paste0("Mouse Embryo E18 (10x, RNA+ATAC)\nATAC UMAP"))
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_mouseembryo_atac-umap.png"),
                plot2, device = "png", width = 7, height = 5, units = "in")

plot3 <- Seurat::DimPlot(mbrain, reduction = "umap.wnn",
                         group.by = "label_Savercat", label = TRUE,
                         cols = col_palette, 
                         repel = TRUE, label.size = 2.5)
plot3 <- plot3 + ggplot2::ggtitle(paste0("Mouse Embryo E18 (10x, RNA+ATAC)\nWNN UMAP"))
plot3 <- plot3 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_mouseembryo_wnn-umap.png"),
                plot3, device = "png", width = 7, height = 5, units = "in")


plot3 <- Seurat::DimPlot(mbrain, reduction = "consensusUMAP",
                         group.by = "label_Savercat", label = TRUE,
                         cols = col_palette, 
                         repel = TRUE, label.size = 2.5)
plot3 <- plot3 + ggplot2::ggtitle(paste0("Mouse Embryo E18 (10x, RNA+ATAC)\nConsensus PCA UMAP"))
plot3 <- plot3 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_mouseembryo_consensusPCA-umap.png"),
                plot3, device = "png", width = 7, height = 5, units = "in")

##########

plot3 <- Seurat::DimPlot(mbrain, reduction = "umap.atac",
                         group.by = "seurat_clusters", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot3 <- plot3 + ggplot2::ggtitle(paste0("Mouse Embryo E18 (10x, RNA+ATAC)\nATAC Seurat clusters"))
plot3 <- plot3 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_mouseembryo_atac-clusters.png"),
                plot3, device = "png", width = 7, height = 5, units = "in")


###########

plot1 <- Seurat::FeaturePlot(mbrain, reduction = "umap.atac",
                             features = "pseudotime")
plot1 <- plot1 + ggplot2::ggtitle(paste0("Mouse Embryo E18 (10x, RNA+ATAC): ATAC UMAP\nSlingshot's pseudotime via ATAC"))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_mouseembryo_atac-pseudotime.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <- Seurat::FeaturePlot(mbrain, reduction = "umap",
                             features = "pseudotime")
plot1 <- plot1 + ggplot2::ggtitle(paste0("Mouse Embryo E18 (10x, RNA+ATAC): RNA UMAP\nSlingshot's pseudotime via ATAC"))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_mouseembryo_rna-pseudotime.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot3 <- Seurat::FeaturePlot(mbrain, reduction = "umap.wnn",
                             features = "pseudotime")
plot3 <- plot3 + ggplot2::ggtitle(paste0("Mouse Embryo E18 (10x, RNA+ATAC): WNN UMAP\nSlingshot's pseudotime via ATAC"))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_mouseembryo_wnn-pseudotime.png"),
                plot3, device = "png", width = 6, height = 5, units = "in")

######

plot1 <- Seurat::DimPlot(mbrain, reduction = "umap.atac",
                         cells.highlight = which(mbrain@meta.data[,"Lineage1"] != 0))
plot1 <- plot1 + Seurat::NoLegend()
plot2 <- Seurat::DimPlot(mbrain, reduction = "umap.atac",
                         cells.highlight = which(mbrain@meta.data[,"Lineage2"] != 0))
plot2 <- plot2 + Seurat::NoLegend()
plot3 <- Seurat::DimPlot(mbrain, reduction = "umap.atac",
                         cells.highlight = which(mbrain@meta.data[,"Lineage3"] != 0))
plot3 <- plot3 + Seurat::NoLegend()
plot_all <- cowplot::plot_grid(plot1, plot2, plot3, nrow = 1, ncol = 3)
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_mouseembryo_atac_lineageHighlighted.png"),
                plot_all, device = "png", width = 10, height = 4, units = "in")


mbrain[["umap.atac"]]@cell.embeddings <- tiltedCCA:::.rotate_matrix(
  source_mat = mbrain[["umap"]]@cell.embeddings,
  target_mat = mbrain[["umap.atac"]]@cell.embeddings
)
plot1 <- Seurat::DimPlot(mbrain, reduction = "umap.atac",
                         group.by = "label_Savercat", 
                         cols = col_palette, pt.size = 1.5)
plot1 <- plot1 + Seurat::NoLegend() 
plot1 <- plot1 + ggplot2::ggtitle("") + ggplot2::ylab("ATAC UMAP 2") + ggplot2::xlab("ATAC UMAP 1")
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_mouseembryo_atac-umap_slides.png"),
                plot1, device = "png", width = 4, height = 4, units = "in")

plot1 <- Seurat::DimPlot(mbrain, reduction = "umap",
                         group.by = "label_Savercat", 
                         cols = col_palette, pt.size = 1.5)
plot1 <- plot1 + Seurat::NoLegend() 
plot1 <- plot1 + ggplot2::ggtitle("") + ggplot2::ylab("RNA UMAP 2") + ggplot2::xlab("RNA UMAP 1")
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_mouseembryo_rna-umap_slides.png"),
                plot1, device = "png", width = 4, height = 4, units = "in")

