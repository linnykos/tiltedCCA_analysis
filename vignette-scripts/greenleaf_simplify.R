rm(list=ls())
library(Seurat)
library(Signac)

load("../../../out/main/10x_greenleaf_preprocessed.RData")

# remove gene activity
greenleaf[["geneActivity"]] <- NULL

# remove RNA
greenleaf[["RNA"]] <- NULL

# remove the embeddings
greenleaf[["pca"]] <- NULL
greenleaf[["lsi"]] <- NULL
greenleaf[["umap.atac"]] <- NULL
greenleaf[["umap.wnn"]] <- NULL
greenleaf[["pcaGeneActivity"]] <- NULL
greenleaf[["umap.geneActivity"]] <- NULL
greenleaf[["consensusPCA"]] <- NULL
greenleaf[["consensusUMAP"]] <- NULL
greenleaf[["umap"]] <- NULL

# remove the graphs
greenleaf@graphs$wknn <- NULL
greenleaf@graphs$wsnn <- NULL
greenleaf@graphs$ATAC_nn <- NULL
greenleaf@graphs$ATAC_snn <- NULL

# remove some metadata
greenleaf@meta.data[,c("SCT.weight", 
                       "ATAC.weight", 
                       "nCount_geneActivity", 
                       "nFeature_geneActivity",
                       "ATAC_snn_res.2", 
                       "seurat_clusters",
                       "pseudotime",
                       "Lineage1",
                       "Lineage2",
                       "Lineage3")] <- NULL
greenleaf@meta.data[,c("percent.mt", 
                       "percent.rb",
                       "S.Score",
                       "G2M.Score",
                       "Phase")] <- NULL
greenleaf@meta.data[,c("nCount_SCT", 
                       "nFeature_SCT")] <- NULL

mentioned_genes <- c("ALDH2", "APOE", "AQP4", "ASCL1", "BHLHE22",
                     "BHLHE40", "C16orf89", "CAV2", "DLX2", "DOK5",
                     "DUSP1", "EOMES", "ETV4", "FOS", "FOXJ1", "GLI3",
                     "HAS2", "HES1", "HES4", "HSPA1A", "HSPA1B",
                     "ID3", "IGFBP7", "JUN", "KIF1A", "LIMCH1",
                     "MBP", "MEF2C", "NEUROD1", "NEUROD2", "NEUROD4",
                     "NEUROD6", "NEUROG1", "NEUROG2", "NFIA", "NFIB", "NFIC",
                     "NHLH1", "NR2F1", "PAX6", "RFX4", "RUNX1",
                     "OLIG1", "OLIG2", "SOX2", "SOX3", "SOX6",
                     "SOX9", "SOX10", "SOX21", "SPARCL1", "SNCB", "TBX",
                     "TNC", "TOP2A", "TRB1", "WNT11")
mentioned_genes <- intersect(mentioned_genes, rownames(greenleaf[["SCT"]]@data))
greenleaf[["SCT"]]@var.features <- unique(c(greenleaf[["SCT"]]@var.features, 
                                            mentioned_genes))

# simplify SCT
new_sct <- Seurat::CreateAssayObject(data = greenleaf[["SCT"]]@data[greenleaf[["SCT"]]@var.features,])
new_sct@var.features <- greenleaf[["SCT"]]@var.features

# simplify ATAC
new_atac <- Seurat::CreateAssayObject(data = greenleaf[["ATAC"]]@data[greenleaf[["ATAC"]]@var.features,])
new_atac@var.features <- greenleaf[["ATAC"]]@var.features

Seurat::DefaultAssay(greenleaf) <- "ATAC"
greenleaf[["SCT"]] <- NULL
greenleaf[["SCT"]] <- new_sct
Seurat::DefaultAssay(greenleaf) <- "SCT"
greenleaf[["ATAC"]] <- NULL
greenleaf[["ATAC"]] <- new_atac

# # remove metadata from ATAC peaks
# greenleaf[["ATAC"]]@fragments[[3]] <- NULL
# greenleaf[["ATAC"]]@fragments[[2]] <- NULL
# greenleaf[["ATAC"]]@fragments[[1]] <- NULL
# greenleaf[["ATAC"]]@annotation <- NULL
# greenleaf[["ATAC"]]@meta.features <- greenleaf[["ATAC"]]@meta.features[,-c(1,2),drop = F]

# throw out other cell types
keep_vec <- rep(TRUE, ncol(greenleaf))
keep_vec[which(greenleaf$celltype %in% c("EC/Peric.", "IN1", "IN2", "IN3", "mGPC/OPC", "SP"))] <- FALSE
greenleaf$keep <- keep_vec
greenleaf <- subset(greenleaf, keep == TRUE)
greenleaf$keep <- NULL

# throw out other cells
seurat_obj <- greenleaf
load("../../../out/main/10x_greenleaf_tcca_RNA-ATAC.RData")

keep_vec <- rep(1, ncol(greenleaf))
keep_vec[greenleaf$celltype %in% c("EC/Peric.", "IN1", "IN2", "IN3", "mGPC/OPC", "SP")] <- 0
greenleaf$keep <- keep_vec
greenleaf <- subset(greenleaf, keep == 1)

greenleaf[["umap.atac"]]@cell.embeddings <- tiltedCCA:::.rotate_matrix(
  source_mat = greenleaf[["umap"]]@cell.embeddings,
  target_mat = greenleaf[["umap.atac"]]@cell.embeddings
)

keep_vec <- rep(1, ncol(greenleaf))
keep_vec[greenleaf[["umap"]]@cell.embeddings[,2]>=5] <- 0
keep_vec[greenleaf[["umap"]]@cell.embeddings[,2]<=-11] <- 0
keep_vec[greenleaf[["umap.atac"]]@cell.embeddings[,2]>=5.5] <- 0
keep_vec[greenleaf[["umap.atac"]]@cell.embeddings[,2]<=-11] <- 0
greenleaf$keep <- keep_vec
greenleaf <- subset(greenleaf, keep == 1)

cell_names <- colnames(greenleaf)
all(cell_names %in% colnames(seurat_obj))
keep_vec <- rep(FALSE, ncol(seurat_obj))
keep_vec[which(colnames(seurat_obj) %in% cell_names)] <- TRUE
seurat_obj$keep <- keep_vec
seurat_obj <- subset(seurat_obj, keep == TRUE)
seurat_obj$keep <- NULL

#########

save(seurat_obj, 
     file = "../../../out/main/10x_trevino_simplified.RData")
