rm(list=ls())

# load in the older version of the Greenleaf data to 1) determine which cells to keep and 2) transfer its cell metadata over
# Technically, we could grab these metadata files from $HOME/nzhanglab/data/GSE162170_multiome_cell_metadata.txt.gz as well. Oops!
load("../../../../out/Writeup14l/Writeup14l_10x_greenleaf_preprocess.RData")

batch_name <- sapply(rownames(greenleaf@meta.data), function(str){
  tmp <- strsplit(str, split = "_")[[1]]
  paste0(tmp[1:(length(tmp)-1)], collapse = "_")
})

library(Seurat); library(Signac)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

create_seurat_object <- function(file_prefix, file_folder, file_suffix){
  tmp <- Seurat::Read10X_h5(paste0(file_prefix, file_folder, file_suffix))
  
  seurat_obj <- Seurat::CreateSeuratObject(
    counts = tmp[["Gene Expression"]],
    assay = "RNA"
  )
  
  seurat_obj
}

##########################

file_prefix <- "~/nzhanglab/data/GSE162170_cortical_multiome/raw-data/"
file_suffix <- "/filtered_feature_bc_matrix.h5"
file_folders <- c("hft_ctx_w21_dc1r3_r1", "hft_ctx_w21_dc2r2_r1", 
                  "hft_ctx_w21_dc2r2_r2")
name_vec <- c("hft_ctx_w21_dc1r3_r1", "hft_ctx_w21_dc2r2_r1", "hft_ctx_w21_dc2r2_r2")

annotation <- Signac::GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
ensembldb::seqlevelsStyle(annotation) <- "UCSC"
Signac::genome(annotation) <- "hg38"

seurat_list <- lapply(1:length(file_folders), function(i){
  print(i)
  file_folder <- file_folders[i]
  seurat_obj <- create_seurat_object(file_prefix = file_prefix,
                                     file_folder = file_folder,
                                     file_suffix = file_suffix)
  seurat_obj$dataset <- name_vec[i]
  print(dim(seurat_obj))
  seurat_obj
})

greenleaf_new <- merge(seurat_list[[1]], y = c(seurat_list[[2]], seurat_list[[3]]), 
                       add.cell.ids = name_vec, 
                       project = "hft", merge.data = T)
table(batch_name)
table(greenleaf_new$dataset)

###################

# read in peak sets
peaks_list <- lapply(file_folders, function(file_folder){
  read.table(
    file = paste0(file_prefix, file_folder, "/atac_peaks.bed"),
    col.names = c("chr", "start", "end")
  )
})

# convert to genomic ranges
gr_list <- lapply(peaks_list, function(peaks_data){
  GenomicRanges::makeGRangesFromDataFrame(peaks_data)
})

# Create a unified set of peaks to quantify in each dataset
combined_peaks <- reduce(x = c(gr_list[[1]], gr_list[[2]], gr_list[[3]]))

# Filter out bad peaks based on length
peak_widths <- IRanges::width(combined_peaks)
combined_peaks <- combined_peaks[peak_widths < 10000 & peak_widths > 20]

# create fragment objects
frags_list <- lapply(1:length(file_folders), function(i){
  file_folder <- file_folders[i]
  Signac::CreateFragmentObject(
    path = paste0(file_prefix, file_folder, "/atac_fragments.tsv.gz"),
    cells = colnames(seurat_list[[i]])
  )
})

atac_count_list <- lapply(1:length(file_folders), function(i){
  Signac::FeatureMatrix(
    fragments = frags_list[[i]],
    features = combined_peaks,
    cells = colnames(seurat_list[[i]])
  )
})

seurat_atac_list <- lapply(1:length(file_folders), function(i){
  assay_obj <- Signac::CreateChromatinAssay(counts = atac_count_list[[i]],
                                            fragments = frags_list[[i]],
                                            sep = c("-", "-"),
                                            annotation = annotation)
  seurat_obj <- Seurat::CreateSeuratObject(assay_obj, assay = "ATAC")
  seurat_obj$dataset <- name_vec[i]
  seurat_obj
})

greenleaf_atac <- merge(seurat_atac_list[[1]], y = c(seurat_atac_list[[2]], seurat_atac_list[[3]]), 
                        add.cell.ids = name_vec, 
                        project = "hft", merge.data = T)

greenleaf_new[["ATAC"]] <- greenleaf_atac[["ATAC"]]
keep_vec <- rep(0, ncol(greenleaf_new))
tmp <- colnames(greenleaf_new); tmp <- sapply(tmp, function(x){strsplit(x, split = "-")[[1]][1]}); names(tmp) <- NULL
keep_vec[which(tmp %in% colnames(greenleaf))] <- 1
table(keep_vec)
greenleaf_new$keep <- keep_vec
greenleaf_new <- subset(greenleaf_new, keep == 1)

########################################3

name_vec <- colnames(greenleaf_new)
name_vec <- sapply(name_vec, function(x){strsplit(x, split = "-")[[1]][1]})
names(name_vec) <- NULL

age_vec <- rep(NA, ncol(greenleaf_new)) # Sample.Age
for(i in 1:ncol(greenleaf_new)){
  age_vec[i] <- greenleaf$Sample.Age[which(rownames(greenleaf@meta.data) == name_vec[i])]
}
table(age_vec, greenleaf$Sample.Age)

celltype_vec <- rep(NA, ncol(greenleaf_new)) # celltype
for(i in 1:ncol(greenleaf_new)){
  celltype_vec[i] <- greenleaf$celltype[which(rownames(greenleaf@meta.data) == name_vec[i])]
}
table(celltype_vec, greenleaf$celltype)

batch_vec <- rep(NA, ncol(greenleaf_new)) # Sample.Batch
for(i in 1:ncol(greenleaf_new)){
  batch_vec[i] <- greenleaf$Sample.Batch[which(rownames(greenleaf@meta.data) == name_vec[i])]
}
table(batch_vec, greenleaf$Sample.Batch)

greenleaf_new$Sample.Age <- age_vec
greenleaf_new$celltype <- celltype_vec
greenleaf_new$Sample.Batch <- batch_vec

###################

set.seed(10)
Seurat::DefaultAssay(greenleaf_new) <- "RNA"
greenleaf_new <- Seurat::SCTransform(greenleaf_new)
greenleaf_new <- Seurat::FindVariableFeatures(greenleaf_new)
greenleaf_new <- Seurat::RunPCA(greenleaf_new, verbose = FALSE)
set.seed(10)
greenleaf_new <- Seurat::RunUMAP(greenleaf_new, dims = 1:50)

set.seed(10)
DefaultAssay(greenleaf_new) <- "ATAC"
greenleaf_new <- Signac::RunTFIDF(greenleaf_new)
greenleaf_new <-  Signac::FindTopFeatures(greenleaf_new, min.cutoff="q10")
greenleaf_new <-  Signac::RunSVD(greenleaf_new)  
set.seed(10)
greenleaf_new <- Seurat::RunUMAP(greenleaf_new, reduction="lsi", 
                             dims=2:50, reduction.name="umap.atac", 
                             reduction.key="atacUMAP_")

####################

set.seed(10)
greenleaf_new <- Seurat::FindMultiModalNeighbors(greenleaf_new, reduction.list = list("pca", "lsi"), 
                                             dims.list = list(1:50, 2:50))
set.seed(10)
greenleaf_new <- Seurat::RunUMAP(greenleaf_new, nn.name = "weighted.nn", reduction.name = "umap.wnn", 
                             reduction.key = "wnnUMAP_")

Seurat::DefaultAssay(greenleaf_new) <- "SCT"
greenleaf_new[["percent.mt"]] <- Seurat::PercentageFeatureSet(object = greenleaf_new, pattern = "^MT-")
greenleaf_new[["percent.rb"]] <- Seurat::PercentageFeatureSet(object = greenleaf_new, pattern = "^RPS")
greenleaf_new <- Seurat::CellCycleScoring(greenleaf_new, 
                                      g2m.features = cc.genes$g2m.genes, 
                                      s.features = cc.genes$s.genes)

###################

plot1 <- Seurat::DimPlot(greenleaf_new, reduction = "umap",
                         group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nRNA UMAP"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14n/Writeup14n_greenleaf_rna-umap.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")


plot2 <- Seurat::DimPlot(greenleaf_new, reduction = "umap.atac",
                         group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot2 <- plot2 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nATAC UMAP"))
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14n/Writeup14n_greenleaf_atac-umap.png"),
                plot2, device = "png", width = 6, height = 5, units = "in")

plot3 <- Seurat::DimPlot(greenleaf_new, reduction = "umap.wnn",
                         group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot3 <- plot3 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nWNN UMAP"))
plot3 <- plot3 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14n/Writeup14n_greenleaf_wnn-umap.png"),
                plot3, device = "png", width = 6, height = 5, units = "in")

####################

# load in gene activity
load("~/nzhanglab/data/GSE162170_cortical_multiome/GSE162170_multiome_atac_gene_activities.RData")

name_vec <- colnames(greenleaf_new)
name_vec <- sapply(name_vec, function(x){strsplit(x, split = "-")[[1]][1]})
names(name_vec) <- NULL

mat <- mat[,name_vec]
colnames(mat) <- colnames(greenleaf_new)
rownames(mat) <- paste0("ATAC-", rownames(mat))

greenleaf_new[["geneActivity"]] <- Seurat::CreateAssayObject(counts = mat)

###################

date_of_run <- Sys.time()
session_info <- devtools::session_info()

greenleaf <- greenleaf_new
save(greenleaf, date_of_run, session_info,
     file = "../../../../out/Writeup14n/Writeup14n_10x_greenleaf.RData")


