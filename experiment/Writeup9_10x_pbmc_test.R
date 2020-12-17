## mainly from: https://satijalab.org/signac/articles/pbmc_multiomic.html
rm(list=ls())

library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(dplyr)
library(ggplot2)

dat <- Matrix::readMM("data/filtered_feature_bc_matrix/matrix.mtx")
dat2 <- as.matrix(dat)

barcodes <- read.table("data/filtered_feature_bc_matrix/barcodes.tsv", sep = "\t")
features <- read.table("data/filtered_feature_bc_matrix/features.tsv", sep = "\t")
head(features)
nrow(features)
table(features[,3])
table(features[,4])

colnames(dat2) <- barcodes[,1]
rownames(dat2) <- features[,1]

# extract RNA and ATAC data
gene_idx <- which(features[,3] == "Gene Expression")
atac_idx <- which(features[,3] == "Peaks")
rna_counts <- dat2[gene_idx,]
atac_counts <- dat2[atac_idx,]

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
frag_file <- "data/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"
chrom_assay <- Signac::CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = frag_file,
  min.cells = 10,
  annotation = annotations
)

pbmc[["ATAC"]] <- chrom_assay
pbmc <- base::subset(
  x = pbmc,
  subset = nCount_ATAC < 7e4 &
    nCount_ATAC > 5e3 &
    nCount_RNA < 25000 &
    nCount_RNA > 1000 &
    percent.mt < 20
)

## perform pre-processing and dimensional reduction on both assays independently
## using standard approaches for rna and atac-seq data
Seurat::DefaultAssay(pbmc) <- "RNA"
pbmc <- Seurat::SCTransform(pbmc, verbose = FALSE)
pbmc <- Seurat::RunPCA(pbmc) 
pbmc <- Seurat::RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

Seurat::DefaultAssay(pbmc) <- "ATAC"
pbmc <- Signac::FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- Signac::RunTFIDF(pbmc)
pbmc <- Signac::RunSVD(pbmc)
pbmc <- Seurat::RunUMAP(pbmc, reduction = 'lsi', dims = 2:50,
                        reduction.name = "umap.atac",
                        reduction.key = "atacUMAP_")

object.size(pbmc)/10^9

tmp <- pbmc[["umap.rna"]]
umap_rna <- cbind(tmp[[,1]], tmp[[,2]])
tmp <- pbmc[["umap.atac"]]
umap_atac <- cbind(tmp[[,1]], tmp[[,2]])

save(umap_rna, umap_atac, file = "out/PBMC_10x_umap.RData")
