rm(list=ls())
library(Seurat)
load("CITEseq_bone_marrow_25antibody.RData")
ls()
bm
head(bm@meta.data) # the metadata, including the cell-types in celltype.l2
table(bm$celltype.l2)

#RNA data
bm[["RNA"]]@counts[1:5,1:5] # the raw count data
dim(bm[["RNA"]]@counts)

head(bm[["RNA"]]@var.features) # the 2000 genes after filtering
length(bm[["RNA"]]@var.features)

bm[["RNA"]]@scale.data[1:5,1:5] # the preprocessed, rescaled/recentered data
dim(bm[["RNA"]]@scale.data)

head(bm[["rna.umap"]]@cell.embeddings) # UMAP coordinates for RNA

##########

#protein data
bm[["ADT"]]@counts[1:5,1:5] # the raw count data
dim(bm[["ADT"]]@counts)

bm[["ADT"]]@scale.data[1:5,1:5] # the preprocessed, rescaled/recentered data
dim(bm[["ADT"]]@scale.data)

head(bm[["adt.umap"]]@cell.embeddings) # UMAP coordinates for ADT

##########

# make some plots
plot1 <- Seurat::DimPlot(bm, reduction = "rna.umap",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human Bone marrow\nCITE-Seq (RNA)"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("CITEseq_bone_marrow_25antibody_rna_umap.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <- Seurat::DimPlot(bm, reduction = "adt.umap",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human Bone marrow\nCITE-Seq (Protein)"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("CITEseq_bone_marrow_25antibody_adt_umap.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")


