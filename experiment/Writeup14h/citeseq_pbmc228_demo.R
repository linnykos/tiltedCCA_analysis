rm(list=ls())
library(Seurat)
load("CITEseq_pbmc_224antibody.RData")
ls()
pbmc
head(pbmc@meta.data) # the metadata, including the cell-types in celltype.l2
table(pbmc$celltype.l2)

#RNA data
pbmc[["SCT"]]@counts[1:5,1:5] # the raw count data
dim(pbmc[["SCT"]]@counts)

head(pbmc[["SCT"]]@var.features) # the 5000 genes after filtering
length(pbmc[["SCT"]]@var.features)

pbmc[["SCT"]]@scale.data[1:5,1:5] # the preprocessed, rescaled/recentered data
dim(pbmc[["SCT"]]@scale.data)

head(pbmc[["rna.umap"]]@cell.embeddings) # UMAP coordinates for RNA

##########

#protein data
pbmc[["ADT"]]@counts[1:5,1:5] # the raw count data
dim(pbmc[["ADT"]]@counts)

pbmc[["ADT"]]@scale.data[1:5,1:5] # the preprocessed, rescaled/recentered data
dim(pbmc[["ADT"]]@scale.data)

head(pbmc[["adt.umap"]]@cell.embeddings) # UMAP coordinates for ADT

##########

# make some plots
plot1 <- Seurat::DimPlot(pbmc, reduction = "rna.umap",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster = F)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human PBMC\nCITE-Seq (RNA)"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("CITEseq_pbmc_224antibody_rna_umap.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <- Seurat::DimPlot(pbmc, reduction = "adt.umap",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster = F)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human PBMC\nCITE-Seq (Protein)"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("CITEseq_pbmc_224antibody_adt_umap.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")


