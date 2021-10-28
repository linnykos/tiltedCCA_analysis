rm(list=ls())
library(Seurat)
library(Signac)

histone_names <- c("H3K27me3")
load(paste0("../../../../data/Pairedtag_mousebrain_RNA-Histone/seurat_",
            histone_names, ".RData"))
mat_1 <- pairedtag[["pca"]]@cell.embeddings
mat_2 <- pairedtag[["lsi"]]@cell.embeddings
celltype <- pairedtag[["celltype"]]
cell_types <- c("Astro_Myoc", "Astro_Nnat", "CA1", "CA23", "CGE",
                "CT", "DG", "Endothelial", "Ependymal", "L23",
                "L4", "L5", "L6", "Microglia", "NP",
                "Oligo_MFOL", "Oligo_MOL", "OPC", "PT", "Pvalb",
                "Sst", "Subiculum")
color_vec <- c("#FF6B2C", "#FF8817", "#5E801F", "#06C65B", "#EA6FE6",
               "#294646", "#00803D", "#B26F13", "#804F0F", "#005D7F",
               "#008FC6", "#05A8EB", "#547180", "#EB8E1A", "#487F80",
               "#EB523A", "#B2402D", "#7F2F23", "#6EC6C6", "#803E7E",
               "#B255B0", "#AAEB32")
color_df <- data.frame(celltype = cell_types, color = color_vec)

save(mat_1, mat_2, celltype, color_df, file = "../../../../out/Writeup14e_simulation/10272021_H3K27me3_original.RData")
