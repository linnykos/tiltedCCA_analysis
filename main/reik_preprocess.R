rm(list=ls())
library(Seurat)
library(Signac)
library(ArchR)

# reik <- readRDS("~/nzhanglab/data/GSE205117_Reik_mouseembryo/ftp1.babraham.ac.uk/data/processed/rna/seurat.rds")
# metadata <- read.csv("~/nzhanglab/data/GSE205117_Reik_mouseembryo/ftp1.babraham.ac.uk/sample_metadata.txt",
#                      sep = "\t")
# 
# set.seed(10)
# date_of_run <- Sys.time()
# session_info <- devtools::session_info()
# 
# n <- ncol(reik)
# mapping_vec <- rep(NA, n)
# for(i in 1:n){
#   if(i %% floor(n/10) == 0) cat('*')
#   tmp <- intersect(which(metadata$barcode == reik$barcode[i]),
#                    which(metadata$sample == reik$sample[i]))
#   
#   if(length(tmp) > 1) stop()
#   if(length(tmp) > 0) mapping_vec[i] <- tmp 
# }
# 
# reik$stage <- metadata$stage[mapping_vec]
# reik$genotype <- metadata$genotype[mapping_vec]
# reik$pass_rnaQC <- metadata$pass_rnaQC[mapping_vec]
# reik$doublet_score <- metadata$doublet_score[mapping_vec]
# reik$doublet_call <- metadata$doublet_call[mapping_vec]
# reik$celltype <- metadata$celltype[mapping_vec]
# reik$celltype.score <- metadata$celltype.score[mapping_vec]
# reik$TSSEnrichment_atac <- metadata$TSSEnrichment_atac[mapping_vec]
# reik$ReadsInTSS_atac <- metadata$ReadsInTSS_atac[mapping_vec]
# reik$PromoterRatio_atac <- metadata$PromoterRatio_atac[mapping_vec]
# reik$NucleosomeRatio_atac <- metadata$NucleosomeRatio_atac[mapping_vec]
# reik$BlacklistRatio_atac <- metadata$BlacklistRatio_atac[mapping_vec]
# reik$pass_atacQC <- metadata$pass_atacQC[mapping_vec]
# reik$celltype.predicted <- metadata$celltype.predicted[mapping_vec]
# 
# set.seed(10)
# Seurat::DefaultAssay(reik) <- "RNA"
# reik <- Seurat::NormalizeData(reik)
# reik <- Seurat::FindVariableFeatures(reik)
# reik <- Seurat::ScaleData(reik) 
# reik <- Seurat::RunPCA(reik, verbose = F)
# set.seed(10)
# reik <- Seurat::RunUMAP(reik, dims = 1:50)
# 
# save(reik, date_of_run, session_info,
#      file = "../../../out/main/10x_reik_preprocessed.RData")

############

load("../../../out/main/10x_reik_preprocessed.RData")
reik_atac <- ArchR::loadArchRProject(path = "~/nzhanglab/data/GSE205117_Reik_mouseembryo/ftp1.babraham.ac.uk/data/processed/atac/archR")
reik_atac <- ArchR::addPeakMatrix(reik_atac, binarize = FALSE, force = TRUE) 
ArchR::getAvailableMatrices(reik_atac) 

# This is a SummarizedExperiment object
se_obj <- ArchR::getMatrixFromProject(ArchRProj = reik_atac,
                                      useMatrix = "PeakMatrix") # this takes about 25 minutes

peak_granges <- se_obj@colData
n <- ncol(reik)
mapping_vec <- rep(NA, n)
for(i in 1:n){
  if(i %% floor(n/10) == 0) cat('*')
  tmp <- intersect(which(peak_granges$barcode == reik$barcode[i]),
                   which(peak_granges$sample == reik$sample[i]))
  
  if(length(tmp) > 1) stop()
  if(length(tmp) > 0) mapping_vec[i] <- tmp 
}

keep_vec <- rep(1, ncol(reik))
keep_vec[which(is.na(mapping_vec))] <- 0
reik$keep <- keep_vec
reik <- subset(reik, keep == 1)

mapping_vec2 <- mapping_vec[which(!is.na(mapping_vec))]
mat <- se_obj@assays@data$PeakMatrix[,mapping_vec2]
chromosome_vec <- as.character(se_obj@rowRanges@seqnames)
range_vec <- as.character(se_obj@rowRanges@ranges)
p <- length(chromosome_vec)
row_vec <- sapply(1:p, function(j){
  paste0(chromosome_vec[j], "_", range_vec[j])
})
rownames(mat) <- row_vec
reik[["ATAC"]] <- Seurat::CreateAssayObject(counts = mat)

Seurat::DefaultAssay(reik) <- "ATAC"
set.seed(10)
reik <- Signac::RunTFIDF(reik)
reik <-  Signac::FindTopFeatures(reik, min.cutoff="q10")
reik <-  Signac::RunSVD(reik)  

set.seed(10)
reik <- Seurat::RunUMAP(reik, reduction="lsi", 
                        dims=2:50, reduction.name="umap.atac", 
                        reduction.key="atacUMAP_")

save(reik, peak_granges, date_of_run, session_info,
     file = "../../../out/main/10x_reik_preprocessed.RData")


