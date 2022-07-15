rm(list=ls())
library(Seurat)
library(Signac)
library(ArchR)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

# reik_atac <- ArchR::loadArchRProject(path = "~/nzhanglab/data/GSE205117_Reik_mouseembryo/ftp1.babraham.ac.uk/data/processed/atac/archR")
# reik_atac <- ArchR::addPeakMatrix(reik_atac, binarize = FALSE, force = TRUE) # I found this from https://github.com/rargelaguet/mouse_organogenesis_10x_multiome_publication/blob/main/atac/archR/peak_calling/filter_peaks_archR.R
# ArchR::getAvailableMatrices(reik_atac) 
# 
# mat <- ArchR::getMatrixFromProject(ArchRProj = reik_atac,
#                                    useMatrix = "PeakMatrix") # this takes about 25 minutes
# save(mat, date_of_run, session_info,
#      readme,
#      file = "~/nzhanglab/data/GSE205117_Reik_mouseembryo/atac_count.RData")

load("~/nzhanglab/data/GSE205117_Reik_mouseembryo/atac_count.RData")
load("../../../../out/Writeup14o/Writeup14o_reik_preprocessed.RData")

###########

dim(mat@assays@data$PeakMatrix)
colnames(mat@assays@data$PeakMatrix)[1:10]
head(mat@colData$barcode)

metadata <- mat@colData
n <- ncol(reik)
mapping_vec <- rep(NA, n)
for(i in 1:n){
  if(i %% floor(n/10) == 0) cat('*')
  tmp <- intersect(which(metadata$barcode == reik$barcode[i]),
                   which(metadata$sample == reik$sample[i]))
  
  if(length(tmp) > 1) stop()
  if(length(tmp) > 0) mapping_vec[i] <- tmp 
}

keep_vec <- rep(1, ncol(reik))
keep_vec[which(is.na(mapping_vec))] <- 0
reik$keep <- keep_vec
reik <- subset(reik, keep == 1)

mapping_vec2 <- mapping_vec[which(!is.na(mapping_vec))]
mat2 <- mat@assays@data$PeakMatrix[,mapping_vec2]
chromosome_vec <- as.character(mat@rowRanges@seqnames)
range_vec <- as.character(mat@rowRanges@ranges)
p <- length(chromosome_vec)
row_vec <- sapply(1:p, function(j){
  paste0(chromosome_vec[j], "_", range_vec[j])
})
rownames(mat2) <- row_vec
reik[["ATAC"]] <- Seurat::CreateAssayObject(counts = mat2)

Seurat::DefaultAssay(reik) <- "ATAC"
set.seed(10)
reik <- Signac::RunTFIDF(reik)
reik <-  Signac::FindTopFeatures(reik, min.cutoff="q10")
reik <-  Signac::RunSVD(reik)  

set.seed(10)
reik <- Seurat::RunUMAP(reik, reduction="lsi", 
                        dims=2:50, reduction.name="umap.atac", 
                        reduction.key="atacUMAP_")

############

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

plot1 <- Seurat::DimPlot(reik, reduction = "umap.atac",
                         group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Mouse embryo (10x, RNA+ATAC)\nATAC UMAP"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14o/Writeup14o_10x_reik_atac-umap_celltype.png"),
                plot1, device = "png", width = 10, height = 5, units = "in")

plot1 <- Seurat::DimPlot(reik, reduction = "umap.atac",
                         group.by = "orig.ident", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Mouse embryo (10x, RNA+ATAC)\nATAC UMAP"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14o/Writeup14o_10x_reik_atac-umap_orig-ident.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")








