## mainly from: https://satijalab.org/signac/articles/pbmc_multiomic.html
rm(list=ls()); set.seed(10); gcinfo(TRUE)

library(Seurat); library(Signac); library(EnsDb.Hsapiens.v86); 
library(BSgenome.Hsapiens.UCSC.hg38); library(GenomeInfoDb)
library(dplyr); library(ggplot2)
library(SeuratDisk) ## https://mojaveazure.github.io/seurat-disk/index.html
library(reticulate)
reticulate::py_module_available("umap")

date_of_run <- Sys.time(); session_info <- sessionInfo()

# load the RNA and ATAC data
counts <- Seurat::Read10X_h5("../../data/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")
fragpath <- "../../data/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"

# get gene annotations for hg38
annotation <- Signac::GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
GenomeInfoDb::seqlevelsStyle(annotation) <- "UCSC"
Signac::genome(annotation) <- "hg38"

# create a Seurat object containing the RNA adata
pbmc <- Seurat::CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# create ATAC assay and add it to the object
pbmc[["ATAC"]] <- Signac::CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)

#########

DefaultAssay(pbmc) <- "ATAC"

pbmc <- Signac::NucleosomeSignal(pbmc)
pbmc <- Signac::TSSEnrichment(pbmc)

# filter out low quality cells
pbmc <- subset(
  x = pbmc,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1000 &
    nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)
pbmc

############################

# call peaks using MACS2
# installation from https://macs3-project.github.io/MACS/docs/INSTALL.html
# roughly 15 minutes
set.seed(10)
peaks <- Signac::CallPeaks(pbmc, macs2.path = "../../venv36/bin/macs3")

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- GenomeInfoDb::keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- IRanges::subsetByOverlaps(x = peaks, ranges = Signac::blacklist_hg38_unified, invert = TRUE)

# quantify counts in each peak
# note: this takes a long time (roughly 15 minutes?)
macs2_counts <- Signac::FeatureMatrix(
  fragments = Signac::Fragments(pbmc),
  features = peaks,
  cells = colnames(pbmc)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
pbmc[["peaks"]] <- Signac::CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)

###################

rm(list = c("annotation", "counts", "fragpath", "macs2_counts", "peaks")); gc(T)

set.seed(10)
Seurat::DefaultAssay(pbmc) <- "RNA"
pbmc <- Seurat::SCTransform(pbmc, assay = "RNA", verbose = T)
pbmc <- Seurat::RunPCA(pbmc, assay = "SCT") 

##################

Seurat::DefaultAssay(pbmc) <- "peaks"
set.seed(10)
pbmc <- Signac::FindTopFeatures(pbmc, min.cutoff = '5')
pbmc <- Signac::RunTFIDF(pbmc)
pbmc <- Signac::RunSVD(pbmc)

##

pbmc <- Seurat::ScaleData(pbmc, features = Seurat::VariableFeatures(object = pbmc))

######################

# load PBMC reference
reference <- SeuratDisk::LoadH5Seurat("../../data/pbmc_multimodal.h5seurat")

DefaultAssay(pbmc) <- "SCT"

set.seed(10)
# transfer cell type labels from reference to query
transfer_anchors <- Seurat::FindTransferAnchors(
  reference = reference,
  query = pbmc,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50
)

set.seed(10)
predictions <- Seurat::TransferData(
  anchorset = transfer_anchors, 
  refdata = reference$celltype.l2,
  weight.reduction = pbmc[['pca']],
  dims = 1:50
)

pbmc <- Seurat::AddMetaData(
  object = pbmc,
  metadata = predictions
)

# set the cell identities to the cell type predictions
Seurat::Idents(pbmc) <- "predicted.id"

# # set a reasonable order for cell types to be displayed when plotting
# levels(pbmc) <- c("CD4 Naive", "CD4 TCM", "CD4 CTL", "CD4 TEM", "CD4 Proliferating",
#                   "CD8 Naive", "dnT",
#                   "CD8 TEM", "CD8 TCM", "MAIT", "NK", "NK_CD56bright",
#                   "NK Proliferating", "gdT",
#                   "Treg", "B naive", "B intermediate", "B memory", "Plasmablast",
#                   "CD14 Mono", "CD16 Mono",
#                   "cDC1", "cDC2", "pDC", "HSPC", "Eryth", "ASDC", "ILC")

#################################

# build a joint neighbor graph using both assays
set.seed(10)
pbmc <- Seurat::FindMultiModalNeighbors(
  object = pbmc,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

# build a joint UMAP visualization
# need to run: reticulate::conda_install("r-reticulate", "umap-learn") 
# i think it's because it's for miniconda?
set.seed(10)
pbmc <- Seurat::RunUMAP(object = pbmc, nn.name = "wknn", verbose = TRUE)

###################

source_code <- readLines("experiment/Writeup10_10x_pbmc_preprocess3.R")
save.image(file = "../../out/Writeup10_10x_pbmc_preprocess3.RData")
