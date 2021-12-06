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
counts <- Seurat::Read10X_h5("../../../../data/10x_mouseembryo/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix.h5")
rownames(counts[[1]]) <- toupper(rownames(counts[[1]]))
mbrain <- Seurat::CreateSeuratObject(counts=counts[["Gene Expression"]])

# Create an assay for the ATAC
# only use peaks in standard chromosomes.
grange.counts <- Signac::StringToGRanges(rownames(counts$Peaks), sep=c(":","-"))
grange.use <- GenomeInfoDb::seqnames(grange.counts) %in% GenomeInfoDb::standardChromosomes(grange.counts)
atac_counts <- counts$Peaks[as.vector(grange.use),]
annotations <- Signac::GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79) # takes a while

# This step is VERY important! Set the annotations levels style so that
# "counts" and "annotations" both use the same name for chromosomes.
# that is, "chr1" and not "1".
GenomeInfoDb::genome(annotations) <- "mm10"
GenomeInfoDb::seqlevelsStyle(annotations)<-"UCSC" # chromosome 1 is called "chr1"
frag.file <- file.path("../../../../data/10x_mouseembryo/e18_mouse_brain_fresh_5k_atac_fragments.tsv.gz")

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
# Get transferred cell labels from Jane.
# --------------------------------------------------

mbrain.jane <- readRDS("../../../../data/10x_mouseembryo/data_tenx_labels_Jane.rds")
mbrain$label_Savercat <- mbrain.jane$savercatLable

# --------------------------------------------------
# Cluster cells on basis of their scRNA-seq profiles
# --------------------------------------------------
set.seed(10)
mbrain <- Seurat::SCTransform(mbrain)

# include housekeeping genes
var_genes <- mbrain[["SCT"]]@var.features
hk_genes <- read.csv("../../../../data/housekeeping.txt")[,1]
s_genes <- cc.genes$s.genes
g2m_genes <- cc.genes$g2m.genes
# 
# all_genes <- rownames(mbrain[["RNA"]]@counts)
# all_genes <- toupper(all_genes)
# any(hk_genes %in% all_genes)
# any(s_genes %in% all_genes)
# any(g2m_genes %in% all_genes)
# hk_genes <- hk_genes[hk_genes %in% rownames(mbrain[["RNA"]]@counts)]
# 
# mbrain <- Seurat::FindVariableFeatures(mbrain)
# mbrain <- Seurat::RunPCA(mbrain, verbose = FALSE)
# 
# # --------------------------------------------------
# # Cluster cells on basis of their scATAC-seq profiles
# # --------------------------------------------------
# set.seed(10)
# DefaultAssay(mbrain)="ATAC"
# mbrain <- Signac::RunTFIDF(mbrain)
# mbrain <-  Signac::FindTopFeatures(mbrain, min.cutoff="q10")
# mbrain <-  Signac::RunSVD(mbrain)  # question: is this svd run with only the top variable features?
# 
# #############################
# 
# metadata <- mbrain@meta.data
# cell_idx <- which(metadata$label_Savercat %in% c("Oligodendrocyte", "Radial glia", 
#                                                  "Forebrain GABAergic", "Neuroblast", 
#                                                  "Glioblast", "Cortical or hippocampal glutamatergic"))
# 
# n <- nrow(metadata)
# keep_vec <- rep(0, n)
# keep_vec[cell_idx] <- 1
# mbrain[["keep"]] <- keep_vec
# mbrain2 <- subset(mbrain, keep == 1)
# 
# Seurat::DefaultAssay(mbrain) <- "SCT"
# set.seed(10)
# mbrain2 <- Seurat::RunUMAP(mbrain2, reduction = 'pca', dims = 1:30, assay = 'SCT', 
#                            reduction.name = 'umap', reduction.key = 'sctUMAP_')
# set.seed(10)
# mbrain2 <- RunUMAP(mbrain2, reduction = 'lsi', dims = 2:50, assay = 'ATAC', 
#                    reduction.name = 'umap.atac', reduction.key = 'atacUMAP_')
# 
# set.seed(10)
# mbrain2 <- Seurat::FindMultiModalNeighbors(
#   mbrain2,
#   reduction.list = list("pca", "lsi"), 
#   dims.list = list(1:30, 2:50))
# set.seed(10)
# mbrain2 <- Seurat::RunUMAP(mbrain2, 
#                            nn.name = "weighted.nn", 
#                            reduction.name = "wnn.umap", 
#                            reduction.key = "wnnUMAP_")
# 
# 

