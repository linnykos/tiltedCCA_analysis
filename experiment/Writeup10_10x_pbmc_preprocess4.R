rm(list=ls()); set.seed(10); gcinfo(TRUE)

library(Seurat); library(Signac); library(EnsDb.Hsapiens.v86); 
library(BSgenome.Hsapiens.UCSC.hg38); library(GenomeInfoDb)
library(dplyr); library(ggplot2)
library(SeuratDisk) ## https://mojaveazure.github.io/seurat-disk/index.html
library(reticulate)
reticulate::py_module_available("umap")

date_of_run <- Sys.time(); session_info <- sessionInfo()

inputdata_10x <- Seurat::Read10X_h5("../../data/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")
rna_counts <- inputdata_10x$`Gene Expression`; atac_counts <- inputdata_10x$Peaks

# Create Seurat object
pbmc <- Seurat::CreateSeuratObject(counts = rna_counts)
pbmc[["percent.mt"]] <- Seurat::PercentageFeatureSet(pbmc, pattern = "^MT-")

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange_counts <- Signac::StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange_use <- seqnames(grange_counts) %in% GenomeInfoDb::standardChromosomes(grange_counts)
atac_counts <- atac_counts[as.vector(grange_use), ]
annotations <- Signac::GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
GenomeInfoDb::seqlevelsStyle(annotations) <- 'UCSC'
Signac::genome(annotations) <- "hg38"

# note: also requires pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi to be in the directory
frag_file <- "../../data/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"
rm(list=c("inputdata_10x", "rna_counts")); gc(verbose = T)
chrom_assay <- Signac::CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = frag_file,
  min.cells = 10,
  annotation = annotations
)
pbmc[["ATAC"]] <- chrom_assay

pbmc <- subset(
  x = pbmc,
  subset = nCount_ATAC < 7e4 &
    nCount_ATAC > 5e3 &
    nCount_RNA < 25000 &
    nCount_RNA > 1000 &
    percent.mt < 20
)

rm(list=c("annotations", "frag_file", "atac_counts", "grange_counts", "grange_use")); gc(verbose = T)

########################

set.seed(10)
Seurat::DefaultAssay(pbmc) <- "RNA"
pbmc <- Seurat::SCTransform(pbmc, assay = "RNA", verbose = T)
pbmc <- Seurat::RunPCA(pbmc, assay = "SCT") 
set.seed(10)
pbmc <- Seurat::RunUMAP(pbmc, dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

#########################

Seurat::DefaultAssay(pbmc) <- "ATAC"
set.seed(10)
pbmc <- Signac::FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- Signac::RunTFIDF(pbmc)
pbmc <- Signac::RunSVD(pbmc)
set.seed(10)
pbmc <- Seurat::RunUMAP(pbmc, reduction = 'lsi', dims = 2:50,
                        reduction.name = "umap.atac",
                        reduction.key = "atacUMAP_")

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

############################

set.seed(10)
pbmc <- Seurat::FindMultiModalNeighbors(pbmc, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
pbmc <- Seurat::RunUMAP(pbmc, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

###############################

rm(list = c("chrom_assay", "reference", "predictions", "transfer_anchors"))
source_code <- readLines("experiment/Writeup10_10x_pbmc_preprocess4.R")
save.image("../../out/Writeup10_10x_pbmc_preprocess4.RData")



