rm(list=ls())
library(Seurat)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

reik <- readRDS("~/nzhanglab/data/GSE205117_Reik_mouseembryo/ftp1.babraham.ac.uk/data/processed/rna/seurat.rds")
metadata <- read.csv("~/nzhanglab/data/GSE205117_Reik_mouseembryo/ftp1.babraham.ac.uk/sample_metadata.txt",
                     sep = "\t")

colnames(reik@meta.data)
colnames(metadata)
head(reik$barcode)
head(metadata$barcode)

n <- ncol(reik)
mapping_vec <- rep(NA, n)
for(i in 1:n){
  if(i %% floor(n/10) == 0) cat('*')
  tmp <- intersect(which(metadata$barcode == reik$barcode[i]),
                   which(metadata$sample == reik$sample[i]))
  
  if(length(tmp) > 1) stop()
  if(length(tmp) > 0) mapping_vec[i] <- tmp 
}

reik$stage <- metadata$stage[mapping_vec]
reik$genotype <- metadata$genotype[mapping_vec]
reik$pass_rnaQC <- metadata$pass_rnaQC[mapping_vec]
reik$doublet_score <- metadata$doublet_score[mapping_vec]
reik$doublet_call <- metadata$doublet_call[mapping_vec]
reik$celltype <- metadata$celltype[mapping_vec]
reik$celltype.score <- metadata$celltype.score[mapping_vec]
reik$TSSEnrichment_atac <- metadata$TSSEnrichment_atac[mapping_vec]
reik$pass_atacQC <- metadata$pass_atacQC[mapping_vec]
reik$celltype.predicted <- metadata$celltype.predicted[mapping_vec]

set.seed(10)
Seurat::DefaultAssay(reik) <- "RNA"
reik <- Seurat::NormalizeData(reik)
reik <- Seurat::FindVariableFeatures(reik)
reik <- Seurat::ScaleData(reik) 
reik <- Seurat::RunPCA(reik, verbose = F)
set.seed(10)
reik <- Seurat::RunUMAP(reik, dims = 1:50)

###############

plot1 <- Seurat::DimPlot(reik, reduction = "umap",
                         group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Mouse embryo (10x, RNA+ATAC)\nRNA UMAP"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14o/Writeup14o_10x_reik_rna-umap_celltype.png"),
                plot1, device = "png", width = 10, height = 5, units = "in")

plot1 <- Seurat::DimPlot(reik, reduction = "umap",
                         group.by = "orig.ident", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Mouse embryo (10x, RNA+ATAC)\nRNA UMAP"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14o/Writeup14o_10x_reik_rna-umap_orig-ident.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

###############

# There are hdf5 files hm...
# The following line doesn't work. It seems like these arrow files only have the fragments...
tmp <- ArchR::getMatrixFromArrow(ArrowFile = "~/nzhanglab/data/GSE205117_Reik_mouseembryo/ftp1.babraham.ac.uk/data/processed/atac/archR/ArrowFiles/E7.5_rep1.arrow",
                                 useMatrix = "PeakMatrix")


reik_atac <- ArchR::loadArchRProject(path = "~/nzhanglab/data/GSE205117_Reik_mouseembryo/ftp1.babraham.ac.uk/data/processed/atac/archR")
reik_atac <- ArchR::addPeakMatrix(reik_atac, binarize = FALSE, force = TRUE) # I found this from https://github.com/rargelaguet/mouse_organogenesis_10x_multiome_publication/blob/main/atac/archR/peak_calling/filter_peaks_archR.R
ArchR::getAvailableMatrices(reik_atac) # from https://www.archrproject.com/articles/Articles/tutorial.html

# https://github.com/GreenleafLab/ArchR/issues/283
mat <- ArchR::getMatrixFromProject(ArchRProj = reik_atac,
                                   useMatrix = "PeakMatrix") # this takes about 25 minutes

readme <- "This is a file process by Kevin on 07/14/2022. mat is a SummarizedExperiment object"
save(mat, date_of_run, session_info,
     readme,
     file = "~/nzhanglab/data/GSE205117_Reik_mouseembryo/atac_count.RData")
