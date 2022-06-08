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
counts <- Seurat::Read10X_h5("~/nzhanglab/data/10x_mouse_embryo/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix.h5")
rownames(counts[["Gene Expression"]]) <- toupper(rownames(counts[["Gene Expression"]]))
mbrain <- Seurat::CreateSeuratObject(counts = counts[["Gene Expression"]])

# Create an assay for the ATAC
# only use peaks in standard chromosomes.
grange.counts <- Signac::StringToGRanges(rownames(counts$Peaks), sep=c(":","-"))
grange.use <- GenomeInfoDb::seqnames(grange.counts) %in% GenomeInfoDb::standardChromosomes(grange.counts)
atac_counts <- counts[["Peaks"]][as.vector(grange.use),]
annotations <- Signac::GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79) # takes a while

# This step is VERY important! Set the annotations levels style so that
# "counts" and "annotations" both use the same name for chromosomes.
# that is, "chr1" and not "1".
GenomeInfoDb::genome(annotations) <- "mm10"
GenomeInfoDb::seqlevelsStyle(annotations)<-"UCSC" # chromosome 1 is called "chr1"
frag.file <- file.path("~/nzhanglab/data/10x_mouse_embryo/e18_mouse_brain_fresh_5k_atac_fragments.tsv.gz")

# needs the e18_mouse_brain_fresh_5k_atac_fragments.tsv.gz.tbi file also
chrom_assay <- Signac::CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)

mbrain[["ATAC"]] <- chrom_assay # add to Seurat object.

# --------------------------------------------------
# Cluster cells on basis of their scRNA-seq profiles
# --------------------------------------------------
set.seed(10)
Seurat::DefaultAssay(mbrain) <- "RNA"
mbrain <- Seurat::SCTransform(mbrain)
mbrain <- Seurat::FindVariableFeatures(mbrain)
mbrain <- Seurat::RunPCA(mbrain, verbose = FALSE)
set.seed(10)
mbrain <- Seurat::RunUMAP(mbrain, dims = 1:50)

# --------------------------------------------------
# Get transferred cell labels from Jane.
# --------------------------------------------------

mbrain.jane <- readRDS("~/project/tiltedCCA/data/10x_mouseembryo/data_tenx_labels_Jane.rds")
mbrain$label_Savercat <- mbrain.jane$savercatLable

# --------------------------------------------------
# Cluster cells on basis of their scATAC-seq profiles
# --------------------------------------------------
set.seed(10)
DefaultAssay(mbrain) <- "ATAC"
mbrain <- Signac::RunTFIDF(mbrain)
mbrain <-  Signac::FindTopFeatures(mbrain, min.cutoff="q10")
mbrain <-  Signac::RunSVD(mbrain)  # question: is this svd run with only the top variable features?
set.seed(10)
mbrain <- Seurat::RunUMAP(mbrain, 
                          reduction="lsi", 
                          dims=2:50,
                          reduction.name="umap.atac", 
                          reduction.key="atacUMAP_")

# --------------------------------------------------
# Calculate a WNN graph
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
mbrain[["percent.mt"]] <- Seurat::PercentageFeatureSet(object = mbrain, pattern = "^MT-")
mbrain[["percent.rb"]] <- Seurat::PercentageFeatureSet(object = mbrain, pattern = "^RPS")
mbrain <- Seurat::CellCycleScoring(mbrain, 
                                   g2m.features = cc.genes$g2m.genes, 
                                   s.features = cc.genes$s.genes)

save(mbrain, date_of_run, session_info,
     file = "../../../out/main/10x_mouseembryo_preprocessed.RData")

####################

plot1 <- Seurat::DimPlot(mbrain, reduction = "umap",
                         group.by = "label_Savercat", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Mouse Embryo E18 (10x, RNA+ATAC)\nRNA UMAP"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_mouseembryo_rna-umap_full.png"),
                plot1, device = "png", width = 8, height = 5, units = "in")

plot2 <- Seurat::DimPlot(mbrain, reduction = "umap.atac",
                         group.by = "label_Savercat", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot2 <- plot2 + ggplot2::ggtitle(paste0("Mouse Embryo E18 (10x, RNA+ATAC)\nATAC UMAP"))
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_mouseembryo_atac-umap_full.png"),
                plot2, device = "png", width = 8, height = 5, units = "in")

plot3 <- Seurat::DimPlot(mbrain, reduction = "umap.wnn",
                         group.by = "label_Savercat", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot3 <- plot3 + ggplot2::ggtitle(paste0("Mouse Embryo E18 (10x, RNA+ATAC)\nWNN UMAP"))
plot3 <- plot3 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_mouseembryo_wnn-umap_full.png"),
                plot3, device = "png", width = 8, height = 5, units = "in")

#########################################

n <- ncol(mbrain)
keep_vec <- rep(0, n)
keep_vec[which(mbrain$label_Savercat %in% c("Oligodendrocyte", "Radial glia", 
                                            "Forebrain GABAergic", "Neuroblast", 
                                            "Glioblast", "Cortical or hippocampal glutamatergic"))] <- 1
mbrain$keep <- keep_vec
mbrain <- subset(mbrain, keep == 1)

DefaultAssay(mbrain) <- "SCT"
mbrain <- Seurat::RunPCA(mbrain, verbose = FALSE)
set.seed(10)
mbrain <- Seurat::RunUMAP(mbrain, dims = 1:50)

set.seed(10)
DefaultAssay(mbrain) <- "ATAC"
mbrain <-  Signac::RunSVD(mbrain)  
set.seed(10)
mbrain <- Seurat::RunUMAP(mbrain, 
                          reduction="lsi", 
                          dims=2:50,
                          reduction.name="umap.atac", 
                          reduction.key="atacUMAP_")

set.seed(10)
mbrain <- Seurat::FindMultiModalNeighbors(mbrain, 
                                          reduction.list = list("pca", "lsi"), 
                                          dims.list = list(1:50, 2:50))
set.seed(10)
mbrain <- Seurat::RunUMAP(mbrain, 
                          nn.name = "weighted.nn", 
                          reduction.name = "umap.wnn", 
                          reduction.key = "wnnUMAP_")

plot1 <- Seurat::DimPlot(mbrain, reduction = "umap",
                         group.by = "label_Savercat", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Mouse Embryo E18 (10x, RNA+ATAC)\nRNA UMAP"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_mouseembryo_rna-umap.png"),
                plot1, device = "png", width = 7, height = 5, units = "in")

plot2 <- Seurat::DimPlot(mbrain, reduction = "umap.atac",
                         group.by = "label_Savercat", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot2 <- plot2 + ggplot2::ggtitle(paste0("Mouse Embryo E18 (10x, RNA+ATAC)\nATAC UMAP"))
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_mouseembryo_atac-umap.png"),
                plot2, device = "png", width = 7, height = 5, units = "in")

plot3 <- Seurat::DimPlot(mbrain, reduction = "umap.wnn",
                         group.by = "label_Savercat", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot3 <- plot3 + ggplot2::ggtitle(paste0("Mouse Embryo E18 (10x, RNA+ATAC)\nWNN UMAP"))
plot3 <- plot3 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_mouseembryo_wnn-umap.png"),
                plot3, device = "png", width = 7, height = 5, units = "in")

##########

Seurat::DefaultAssay(mbrain) <- "ATAC"
set.seed(10)
mbrain <- Seurat::FindNeighbors(mbrain, dims = 2:50, reduction = "lsi")
set.seed(10)
mbrain <- Seurat::FindClusters(mbrain, resolution = 2)

plot3 <- Seurat::DimPlot(mbrain, reduction = "umap.atac",
                         group.by = "seurat_clusters", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot3 <- plot3 + ggplot2::ggtitle(paste0("Mouse Embryo E18 (10x, RNA+ATAC)\nATAC Seurat clusters"))
plot3 <- plot3 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_mouseembryo_atac-clusters.png"),
                plot3, device = "png", width = 7, height = 5, units = "in")

# remove unwanted clusters
keep_vec <- rep(1, ncol(mbrain))
keep_vec[mbrain$seurat_clusters %in% c("5","15","19","6","18","13","14","16", "4")] <- 0
mbrain$keep <- keep_vec
mbrain_tmp <- subset(mbrain, keep == 1)
table(mbrain_tmp$seurat_clusters)

# ad-hoc merge some clusters for Slingshot to work better
mbrain_tmp$seurat_clusters[mbrain_tmp$seurat_clusters == "1"] <- "8"
mbrain_tmp$seurat_clusters[mbrain_tmp$seurat_clusters == "2"] <- "8"
mbrain_tmp$seurat_clusters[mbrain_tmp$seurat_clusters == "20"] <- "9"
mbrain_tmp$seurat_clusters[mbrain_tmp$seurat_clusters == "11"] <- "0"

set.seed(10)
mbrain_tmp <-  Signac::RunSVD(mbrain_tmp)  
dimmat_x <- mbrain_tmp[["lsi"]]@cell.embeddings[,2:50]

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

# extract pseudotime
pseudotime_mat <- slingshot_curve@assays@data$pseudotime
for(i in 1:3){
  # remove pseudotimes for cells not in the correct lineage
  idx <- mbrain_tmp$seurat_clusters
  pseudotime_mat[which(!mbrain_tmp$seurat_clusters %in% slingshot_lin@metadata$lineages[[i]]),i] <- NA
  idx <- which(!is.na(pseudotime_mat[,i]))
  pseudotime_mat[idx,i] <- rank(pseudotime_mat[idx,i])
}
# compute a weighted rank
lineage_length <- sapply(slingshot_lin@metadata$lineages, function(lin){
  length(which(mbrain_tmp$seurat_clusters %in% lin))
})
pseudotime_vec <- apply(pseudotime_mat, 1, function(x){
  vec <- sapply(1:3, function(i){x[i]/lineage_length[i]})
  mean(vec, na.rm = T)
})
pseudotime_vec_full <- rep(NA, ncol(mbrain))
names(pseudotime_vec_full) <- colnames(mbrain)
pseudotime_vec_full[names(pseudotime_vec)] <- pseudotime_vec
mbrain$pseudotime <- pseudotime_vec_full

# add metadata
for(i in 1:3){
  traj_vec <- rep(0, ncol(mbrain))
  names(traj_vec) <- colnames(mbrain)
  idx <- which(!is.na(pseudotime_mat[,i]))
  traj_vec[names(pseudotime_mat[idx,i])] <- 1
  mbrain <- Seurat::AddMetaData(mbrain, metadata = traj_vec, col.name = paste0("Lineage", i))
}

###############

# consensus pca
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
                                           dims_1 = 1:50, dims_2 = 1:50,
                                           dims_consensus = 1:50,
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

##########



save(mbrain, date_of_run, session_info,
     file = "../../../out/main/10x_mouseembryo_preprocessed.RData")

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

plot1 <- Seurat::DimPlot(mbrain, reduction = "consensusUMAP",
                         group.by = "label_Savercat", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Mouse Embryo E18 (10x, RNA+ATAC)\nConsensus PCA's UMAP"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_mouseembryo_consensusPCA-umap.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

