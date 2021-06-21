rm(list=ls())

# following code adapted from Chi-Yun
library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(dplyr)
library(ggplot2)
library(hdf5r)
library(Matrix)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

rna_counts <- readRDS("../../../../data/ICB_mouse/gmat_RNA_counts_ALL_dim35.rds")
atac_counts <- readRDS("../../../../data/ICB_mouse/ALL_processed_FM.rds")
cell_meta <- readRDS("../../../../data/ICB_mouse/metadat.rds")

rna_counts <- rna_counts[,match(rownames(cell_meta), colnames(rna_counts))]
atac_counts <- atac_counts[,match(rownames(cell_meta), colnames(atac_counts))]

set.seed(10)
ind_2000_each=c()
for(ii in unique(cell_meta$Sample)){
  ind_2000_each=c(ind_2000_each, sample(which(cell_meta$Sample==ii), 2000, replace=F))
}

ind_2000_each=sort(ind_2000_each)
rna_counts=rna_counts[, ind_2000_each]
atac_counts=atac_counts[, ind_2000_each]

############

myeloid <- Seurat::CreateSeuratObject(counts = rna_counts)
myeloid[["percent.mt"]] <- Seurat::PercentageFeatureSet(myeloid, pattern = "^mt-")

grange.counts <- Signac::StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- GenomeInfoDb::seqnames(grange.counts) %in% GenomeInfoDb::standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- Signac::GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79)

GenomeInfoDb::seqlevelsStyle(annotations) <- 'UCSC'
GenomeInfoDb::genome(annotations) <- "mm10"
chrom_assay <- Signac::CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'mm10',
  min.cells = 10,
  annotation = annotations
)
myeloid[["ATAC"]] <- chrom_assay

# assign cell-type labels
myeloid$celltype <- myeloid$orig.ident
myeloid$sample <- cell_meta$Sample[ind_2000_each]

# RNA analysis
Seurat::DefaultAssay(myeloid) <- "RNA"
set.seed(10)
myeloid <- Seurat::SCTransform(myeloid)
myeloid <- Seurat::FindVariableFeatures(myeloid)
myeloid <- Seurat::RunPCA(myeloid, verbose = FALSE)
set.seed(10)
myeloid <- Seurat::RunUMAP(myeloid, dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
Seurat::DefaultAssay(myeloid) <- "ATAC"
set.seed(10)
myeloid <- Signac::RunTFIDF(myeloid)
myeloid <- Signac::FindTopFeatures(myeloid, min.cutoff = 'q0')
myeloid <- Signac::RunSVD(myeloid)
set.seed(10)
myeloid <- Seurat::RunUMAP(myeloid, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

# multimodal analysis
set.seed(10)
myeloid <- Seurat::FindMultiModalNeighbors(myeloid, reduction.list = list("pca", "lsi"), 
                                          dims.list = list(1:50, 2:50))
set.seed(10)
myeloid <- Seurat::RunUMAP(myeloid, nn.name = "weighted.nn", reduction.name = "umap.wnn", 
                          reduction.key = "wnnUMAP_")

####################

png("../../../../out/figures/Writeup14/Writeup14_mouseicb_quality.png", 
    height = 1500, width = 2000, units = "px", res = 300)
Seurat::VlnPlot(myeloid, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
        log = TRUE, pt.size = 0) + NoLegend()
graphics.off()

png("../../../../out/figures/Writeup14/Writeup14_mouseicb_seurat_umap_celltype.png", 
    height = 1500, width = 4500, units = "px", res = 300)
p1 <- Seurat::DimPlot(myeloid, reduction = "umap.rna", group.by = "celltype", 
                      label = TRUE, label.size = 2.5, repel = TRUE) + ggplot2::ggtitle("RNA")
p2 <- Seurat::DimPlot(myeloid, reduction = "umap.atac", group.by = "celltype", 
                      label = TRUE, label.size = 2.5, repel = TRUE) + ggplot2::ggtitle("ATAC")
p3 <- Seurat::DimPlot(myeloid, reduction = "umap.wnn", group.by = "celltype", 
                      label = TRUE, label.size = 2.5, repel = TRUE) + ggplot2::ggtitle("WNN")
p1 + p2 + p3 & Seurat::NoLegend() & ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
graphics.off()

png("../../../../out/figures/Writeup14/Writeup14_mouseicb_seurat_umap_sample.png", 
    height = 1500, width = 4500, units = "px", res = 300)
p1 <- Seurat::DimPlot(myeloid, reduction = "umap.rna", group.by = "sample", 
                      label = TRUE, label.size = 2.5, repel = TRUE) + ggplot2::ggtitle("RNA")
p2 <- Seurat::DimPlot(myeloid, reduction = "umap.atac", group.by = "sample", 
                      label = TRUE, label.size = 2.5, repel = TRUE) + ggplot2::ggtitle("ATAC")
p3 <- Seurat::DimPlot(myeloid, reduction = "umap.wnn", group.by = "sample", 
                      label = TRUE, label.size = 2.5, repel = TRUE) + ggplot2::ggtitle("WNN")
p1 + p2 + p3 & Seurat::NoLegend() & ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
graphics.off()

save(myeloid, file = "../../../../out/Writeup14/Writeup14_mouseicb_preprocess.RData")

################

# compare to Chi-Yun's preprocessing for sanity check
rm(list=ls())
library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(dplyr)
library(ggplot2)
library(hdf5r)
library(Matrix)

myeloid <- readRDS("../../../../data/ICB_mouse/Object_seurat_2000_each.rds")
png("../../../../out/figures/Writeup14/Writeup14_mouseicb_seurat_umap_celltype_chiyun.png", 
    height = 1500, width = 3000, units = "px", res = 300)
p1 <- Seurat::DimPlot(myeloid, reduction = "umap.rna", group.by = "celltype", 
                      label = TRUE, label.size = 2.5, repel = TRUE) + ggplot2::ggtitle("RNA")
p2 <- Seurat::DimPlot(myeloid, reduction = "umap.atac", group.by = "celltype", 
                      label = TRUE, label.size = 2.5, repel = TRUE) + ggplot2::ggtitle("ATAC")
p1 + p2 & Seurat::NoLegend() & ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
graphics.off()

