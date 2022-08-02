rm(list=ls())
library(Seurat); library(Signac); library(EnsDb.Hsapiens.v86)

load("../../../out/main/10x_greenleaf_preprocessed.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

greenleaf[["geneActivity"]] <- NULL
greenleaf[["umap.geneActivity"]] <- NULL

set.seed(10)
gene_activities <- Signac::GeneActivity(greenleaf,
                                        extend.upstream = 500,
                                        extend.downstream = 500,
                                        biotypes = "protein_coding")

set.seed(10)
rownames(gene_activities) <- paste0("ATAC-", rownames(gene_activities))
var_features <- greenleaf[["SCT"]]@var.features

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
var_features <- unique(c(var_features, mentioned_genes))

var_features <- intersect(paste0("ATAC-",var_features), rownames(gene_activities))
greenleaf[["customGAct"]] <- Seurat::CreateAssayObject(counts = gene_activities)
Seurat::DefaultAssay(greenleaf) <- "customGAct"
greenleaf <- Seurat::NormalizeData(greenleaf)
greenleaf[["customGAct"]]@var.features <- var_features
greenleaf <- Seurat::ScaleData(greenleaf)
set.seed(10)
greenleaf <- Seurat::RunPCA(greenleaf, verbose = FALSE,
                            reduction.name = "pcaCustomGAct",
                            reduction.key = "PCACustomGAct_")
set.seed(10)
greenleaf <- Seurat::RunUMAP(greenleaf, reduction="pcaCustomGAct", 
                             dims=2:50, reduction.name="umap.customGAct", 
                             reduction.key = "customGActUMAP_")


plot1 <- Seurat::DimPlot(greenleaf, reduction = "umap.customGAct",
                         group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nCust. gene Activity UMAP"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_greenleaf_customGAct-umap.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

save(greenleaf, date_of_run, session_info,
     file = "../../../out/main/10x_greenleaf_preprocessed_customGAct.RData")

