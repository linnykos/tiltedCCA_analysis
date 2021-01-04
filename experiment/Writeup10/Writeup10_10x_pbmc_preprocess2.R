## mainly from: https://satijalab.org/signac/articles/pbmc_multiomic.html
rm(list=ls()); set.seed(10)

library(Seurat); library(Signac); library(EnsDb.Hsapiens.v86); library(GenomeInfoDb)
library(dplyr); library(ggplot2)

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
pbmc <- Seurat::SCTransform(pbmc, assay = "RNA", verbose = FALSE)
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

###########################

# apply clustering from the vignette
# https://satijalab.org/seurat/v4.0/weighted_nearest_neighbor_analysis.html

set.seed(10)
pbmc <- Seurat::FindMultiModalNeighbors(pbmc, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
pbmc <- Seurat::RunUMAP(pbmc, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
pbmc <- Seurat::FindClusters(pbmc, graph.name = "wsnn", algorithm = 3, verbose = FALSE) ## see http://www.ludowaltman.nl/slm/ for SLM (smart local moving)
pbmc <- Seurat::FindSubCluster(pbmc, cluster = 6, graph.name = "wsnn", algorithm = 3)
Seurat::Idents(pbmc) <- "sub.cluster"

pbmc <- Seurat::RenameIdents(pbmc, '19' = 'pDC','20' = 'HSPC','15' = 'cDC')
pbmc <- Seurat::RenameIdents(pbmc, '0' = 'CD14 Mono', '8' ='CD14 Mono', '5' = 'CD16 Mono')
pbmc <- Seurat::RenameIdents(pbmc, '17' = 'Naive B', '11' = 'Intermediate B', '10' = 'Memory B', '21' = 'Plasma')
pbmc <- Seurat::RenameIdents(pbmc, '7' = 'NK')
pbmc <- Seurat::RenameIdents(pbmc, '4' = 'CD4 TCM', '13'= "CD4 TEM", '3' = "CD4 TCM", '16' ="Treg", '1' ="CD4 Naive", '14' = "CD4 Naive")
pbmc <- Seurat::RenameIdents(pbmc, '2' = 'CD8 Naive', '9'= "CD8 Naive", '12' = 'CD8 TEM_1', '6_0' = 'CD8 TEM_2', '6_1' ='CD8 TEM_2')
pbmc <- Seurat::RenameIdents(pbmc, '18' = 'MAIT')
pbmc <- Seurat::RenameIdents(pbmc, '6_2' ='gdT', '6_3' = 'gdT')
pbmc$celltype <- Seurat::Idents(pbmc)

###################

source_code <- readLines("experiment/Writeup10_10x_pbmc_preprocess2.R")
save.image(file = "../../out/Writeup10_10x_pbmc_preprocess.RData")
