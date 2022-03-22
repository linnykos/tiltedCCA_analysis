rm(list=ls())
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(dplyr)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

obj1 <- readRDS("~/nzhanglab/data/GSE162170_cortical_multiome/GSE162170_multiome_rna_counts.rds")
obj1 <- Matrix::Matrix(as.matrix(obj1), sparse = T)
load("~/nzhanglab/data/GSE162170_cortical_multiome/GSE162170_multiome_atac_gene_activities.RData")
obj2 <- mat
all(ncol(obj1) == ncol(obj2))
all(colnames(obj1) == colnames(obj2))

edb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
gene_name <- AnnotationDbi::mapIds(edb,
                                   keys = rownames(obj1),
                                   keytype = "GENEID", column = "GENENAME", multiVals = "first")
names(gene_name) <- NULL
na_idx <- which(is.na(gene_name))
if(length(na_idx) > 0){gene_name <- gene_name[-na_idx]; obj1 <- obj1[-na_idx,,drop=F]}
rownames(obj1) <- gene_name

length(intersect(rownames(obj1), rownames(obj2)))
length(unique(c(rownames(obj1), rownames(obj2))))

gene_names <- sort(intersect(rownames(obj1), rownames(obj2)))
obj1 <- obj1[gene_names,]
obj2 <- obj2[gene_names,]
obj1 <- obj1[order(rownames(obj1)),]
obj2 <- obj2[order(rownames(obj2)),]
all(rownames(obj1) == rownames(obj2))
rownames(obj2) <- paste0("ATAC_", rownames(obj2))

greenleaf <- Seurat::CreateSeuratObject(counts = obj1)
greenleaf[["ATAC"]] <- Seurat::CreateAssayObject(counts = obj2)

set.seed(10)
Seurat::DefaultAssay(greenleaf) <- "RNA"
greenleaf <- Seurat::SCTransform(greenleaf)
greenleaf <- Seurat::FindVariableFeatures(greenleaf)
greenleaf <- Seurat::RunPCA(greenleaf, verbose = FALSE)
set.seed(10)
greenleaf <- Seurat::RunUMAP(greenleaf, dims = 1:50)

####################

set.seed(10)
DefaultAssay(greenleaf) <- "ATAC"
greenleaf <- Seurat::NormalizeData(greenleaf)
greenleaf[["ATAC"]]@var.features <- paste0("ATAC-", greenleaf[["SCT"]]@var.features)
greenleaf <- Seurat::ScaleData(greenleaf)
greenleaf <- Seurat::RunPCA(greenleaf, verbose = FALSE,
                            reduction.name = "lsi",
                            reduction.key = "LSI_")
set.seed(10)
greenleaf <- Seurat::RunUMAP(greenleaf, reduction="lsi", 
                             dims=2:50, reduction.name="umap.atac", 
                             reduction.key="atacUMAP_")

########################################

cell_metadata <- readr::read_delim("~/nzhanglab/data/GSE162170_cortical_multiome/GSE162170_multiome_cell_metadata.txt.gz")
cluster_metadata <- readr::read_delim("~/nzhanglab/data/GSE162170_cortical_multiome/GSE162170_multiome_cluster_names.txt.gz")
cluster_metadata <- as.data.frame(cluster_metadata)
cluster_metadata <- cluster_metadata[which(cluster_metadata$Assay == "Multiome RNA"),]

all(cell_metadata$Cell.ID == colnames(greenleaf))
celltype_vec <- as.character(cell_metadata$seurat_clusters)
for(i in 1:nrow(cluster_metadata)){
  idx <- which(celltype_vec == cluster_metadata[i,"Cluster.ID"])
  celltype_vec[idx] <- cluster_metadata[i,"Cluster.Name"]
}
greenleaf$celltype <- celltype_vec
greenleaf$Sample.Age <- dplyr::pull(cell_metadata, "Sample.Age")
greenleaf$Sample.Batch <- dplyr::pull(cell_metadata, "Sample.Batch")

#####################

set.seed(10)
greenleaf <- Seurat::FindMultiModalNeighbors(greenleaf, reduction.list = list("pca", "lsi"), 
                                             dims.list = list(1:50, 2:50))
set.seed(10)
greenleaf <- Seurat::RunUMAP(greenleaf, nn.name = "weighted.nn", reduction.name = "umap.wnn", 
                             reduction.key = "wnnUMAP_")

Seurat::DefaultAssay(greenleaf) <- "SCT"
greenleaf[["percent.mt"]] <- Seurat::PercentageFeatureSet(object = greenleaf, pattern = "^MT-")
greenleaf[["percent.rb"]] <- Seurat::PercentageFeatureSet(object = greenleaf, pattern = "^RPS")
greenleaf <- Seurat::CellCycleScoring(greenleaf, 
                                      g2m.features = cc.genes$g2m.genes, 
                                      s.features = cc.genes$s.genes)

#####################

plot1 <- Seurat::DimPlot(greenleaf, reduction = "umap",
                         group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+Gene Activity)\nRNA UMAP")) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot2 <- Seurat::DimPlot(greenleaf, reduction = "umap.atac",
                         group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot2 <- plot2 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+Gene Activity)\nGene activity UMAP")) + Seurat::NoLegend()
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot3 <- Seurat::DimPlot(greenleaf, reduction = "umap.wnn",
                         group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot3 <- plot3 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+Gene Activity)\nWNN UMAP"))
plot3 <- plot3 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

plot4 <- cowplot::plot_grid(plot1, plot2, plot3, ncol = 3)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14m/Writeup14m_10x_greenleaf_umap_rnaAndgeneactivity.png"),
                plot4, device = "png", width = 15, height = 5, units = "in")

save(greenleaf, date_of_run, session_info,
     file = "../../../../out/Writeup14m/Writeup14m_10x_greenleaf_preprocess_geneactivity.RData")

