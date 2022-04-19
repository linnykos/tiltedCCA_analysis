rm(list=ls())

library(Seurat); library(Signac); library(EnsDb.Hsapiens.v86); 
library(BSgenome.Hsapiens.UCSC.hg38); library(GenomeInfoDb)
library(dplyr); library(ggplot2)
library(SeuratDisk) ## https://mojaveazure.github.io/seurat-disk/index.html

set.seed(10); gcinfo(TRUE)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

## see directions from https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis.html

inputdata_10x <- Seurat::Read10X_h5("~/project/tiltedCCA/data/10x_PBMC_RNA-ATAC/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")
rna_counts <- inputdata_10x[["Gene Expression"]]
atac_counts <- inputdata_10x[["Peaks"]]

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
frag_file <- "~/project/tiltedCCA/data/10x_PBMC_RNA-ATAC/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"
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
pbmc <- Seurat::RunPCA(pbmc, assay = "SCT",
                       verbose = F) 
set.seed(10)
pbmc <- Seurat::RunUMAP(pbmc, 
                        dims = 1:50, 
                        reduction.name = 'umap.rna', 
                        reduction.key = 'rnaUMAP_')

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

# see directions on https://satijalab.org/seurat/articles/multimodal_reference_mapping.html
# load PBMC reference
reference <- SeuratDisk::LoadH5Seurat("~/project/tiltedCCA/data/10x_PBMC_RNA-ATAC/pbmc_multimodal.h5seurat")

Seurat::DefaultAssay(reference) <- "SCT"
Seurat::DefaultAssay(pbmc) <- "SCT"
set.seed(10)
# transfer cell type labels from reference to query
# see https://github.com/satijalab/seurat/issues/3989
transfer_anchors <- Seurat::FindTransferAnchors(
  reference = reference,
  query = pbmc,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50,
  recompute.residuals = FALSE
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
pbmc <- Seurat::FindMultiModalNeighbors(pbmc, 
                                        reduction.list = list("pca", "lsi"), 
                                        dims.list = list(1:50, 2:50))
set.seed(10)
pbmc <- Seurat::RunUMAP(pbmc, nn.name = "weighted.nn", 
                        reduction.name = "wnn.umap", 
                        reduction.key = "wnnUMAP_")

###############################

save(pbmc, date_of_run, session_info,
     file = "../../../out/main/10x_pbmc_preprocessed.RData")

###################

plot1 <- Seurat::DimPlot(pbmc, reduction = "umap.rna",
                         group.by = "predicted.id", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster=FALSE)
plot1 <- plot1 + ggplot2::ggtitle(paste0("PBMC (10x, RNA+ATAC)\nRNA UMAP"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_pbmc_rna-umap.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot2 <- Seurat::DimPlot(pbmc, reduction = "umap.atac",
                         group.by = "predicted.id", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster=FALSE)
plot2 <- plot2 + ggplot2::ggtitle(paste0("PBMC (10x, RNA+ATAC)\nATAC UMAP"))
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_pbmc_atac-umap.png"),
                plot2, device = "png", width = 6, height = 5, units = "in")

plot3 <- Seurat::DimPlot(pbmc, reduction = "wnn.umap",
                         group.by = "predicted.id", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster=FALSE)
plot3 <- plot3 + ggplot2::ggtitle(paste0("PBMC (10x, RNA+ATAC)\nWNN UMAP"))
plot3 <- plot3 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_pbmc_wnn-umap.png"),
                plot3, device = "png", width = 6, height = 5, units = "in")


