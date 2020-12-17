## mainly from: https://satijalab.org/signac/articles/pbmc_multiomic.html
rm(list=ls())
set.seed(10)

library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(dplyr)
library(ggplot2)

date_of_run <- Sys.time()
session_info <- sessionInfo()

dat <- Matrix::readMM("../../data/filtered_feature_bc_matrix/matrix.mtx")
dat2 <- as.matrix(dat)

barcodes <- read.table("../../data/filtered_feature_bc_matrix/barcodes.tsv", sep = "\t")
features <- read.table("../../data/filtered_feature_bc_matrix/features.tsv", sep = "\t")
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
pbmc@assays$RNA@counts[80:90,90:100]
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

# compare the SCTransform-UMAP and log-normalize-PCA embedding
dim(pbmc@assays$RNA@counts)
Seurat::DefaultAssay(pbmc) <- "RNA"
pbmc2 <- pbmc
pbmc2@assays$RNA@counts[80:90,90:100]
pbmc <- Seurat::SCTransform(pbmc, assay = "RNA", verbose = FALSE)
pbmc@assays$RNA@counts[80:90,90:100]
pbmc <- Seurat::RunPCA(pbmc, assay = "SCT") 
head(pbmc2@reductions$pca)[1:5,1:5]
pbmc <- Seurat::RunUMAP(pbmc, dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

pca_sct_rna <- pbmc[["pca"]]@cell.embeddings[,1:20]
umap_sct_rna <- pbmc[["umap.rna"]]@cell.embeddings

###

Seurat::DefaultAssay(pbmc2) <- "RNA"
pbmc2 <- Seurat::NormalizeData(pbmc2, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc2 <- Seurat::FindVariableFeatures(pbmc2, selection.method = "vst", nfeatures = 2000)
pbmc2 <- Seurat::ScaleData(pbmc2, features = rownames(pbmc2))
pbmc2 <- Seurat::RunPCA(pbmc2, features = Seurat::VariableFeatures(object = pbmc2)) 

tmp <- pbmc2[["RNA"]]@scale.data
tmp[70:80,90:100]
mean(tmp[70,]); sd(tmp[70,]); length(which(abs(tmp[70,]) == 10))
zz <- which(apply(tmp, 1, function(x){any(x != 0) & length(which(abs(x)==10)) == 0}))
head(zz)
mean(tmp[7,]); sd(tmp[7,]); length(which(abs(tmp[7,]) == 10))
pca_log_rna <- pbmc2[["pca"]]@cell.embeddings[,1:20]

#####################

# from https://satijalab.org/signac/articles/pbmc_vignette.html
Seurat::DefaultAssay(pbmc) <- "ATAC"
pbmc <- Signac::FindTopFeatures(pbmc, min.cutoff = 'q0') # there's there's contention on the order of operations here
pbmc <- Signac::RunTFIDF(pbmc)
pbmc <- Signac::RunSVD(pbmc)
pbmc <- Seurat::RunUMAP(pbmc, reduction = 'lsi', dims = 2:50,
                        reduction.name = "umap.atac",
                        reduction.key = "atacUMAP_")

umap_tfidf_atac <- pbmc[["umap.atac"]]@cell.embeddings

### now manually compute the PCA since it's too big to use the usual Seurat way of doing so
pbmc <- Seurat::ScaleData(pbmc, assay = "ATAC", features = rownames(pbmc))
mat <- pbmc[["ATAC"]]@scale.data
mat[70:80,90:100]
mean(mat[1,]); sd(mat[1,]); length(which(abs(mat[1,]) == 10))
zz <- which(apply(mat, 1, function(x){any(x != 0) & length(which(abs(x)==10)) == 0}))
head(zz)
mean(mat[8,]); sd(mat[8,]); length(which(abs(mat[8,]) == 10))
class(mat)
svd_res <- RSpectra::svds(mat, k = 20)

pca_tfidf_atac <- t(mat) %*% svd_res$u 

save(pca_sct_rna, umap_sct_rna, pca_log_rna, umap_tfidf_atac, pca_tfidf_atac,
     file = "../../out/Writeup9_10x_pbmc_dimred_onlyembedding.RData")
rm(list = c("dat", "dat2", "mat", "svd_res", "tmp", "zz", "pbmc2",
            "gene_idx", "atac_idx", "rna_counts", "atac_counts", "chrom_assay",
            "frag_file"))
source_code <- readLines("Writeup9_10x_pbmc_dimred.R")
save.image(file = "../../out/Writeup9_10x_pbmc_dimred_all.RData")
