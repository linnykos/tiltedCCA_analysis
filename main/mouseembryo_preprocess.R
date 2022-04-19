rm(list=ls())

library(Signac)
library(Seurat)
library(EnsDb.Mmusculus.v79)
library(GenomeInfoDb)
library(ggplot2)
library(patchwork)
library(biovizBase)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

# --------------------------------------------------
# Load in data and set up Seurat object
# --------------------------------------------------

# Create a Seurat object first with the RNA
counts <- Seurat::Read10X_h5("~/nzhanglab/data/10x_mouse_embryo/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix.h5")
rownames(counts[["Gene Expression"]]) <- toupper(rownames(counts[["Gene Expression"]]))
mbrain <- Seurat::CreateSeuratObject(counts = counts[["Gene Expression"]])

# Create an assay for the ATAC
# only use peaks in standard chromosomes.
grange.counts <- Signac::StringToGRanges(rownames(counts$Peaks), sep=c(":","-"))
grange.use <- GenomeInfoDb::seqnames(grange.counts) %in% GenomeInfoDb::standardChromosomes(grange.counts)
atac_counts <- counts[["Peaks"]][as.vector(grange.use),]
annotations <- Signac::GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79) # takes a while

# This step is VERY important! Set the annotations levels style so that
# "counts" and "annotations" both use the same name for chromosomes.
# that is, "chr1" and not "1".
GenomeInfoDb::genome(annotations) <- "mm10"
GenomeInfoDb::seqlevelsStyle(annotations)<-"UCSC" # chromosome 1 is called "chr1"
frag.file <- file.path("~/nzhanglab/data/10x_mouse_embryo/e18_mouse_brain_fresh_5k_atac_fragments.tsv.gz")

# needs the e18_mouse_brain_fresh_5k_atac_fragments.tsv.gz.tbi file also
chrom_assay <- Signac::CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)

mbrain[["ATAC"]] <- chrom_assay # add to Seurat object.

# --------------------------------------------------
# Cluster cells on basis of their scRNA-seq profiles
# --------------------------------------------------
set.seed(10)
Seurat::DefaultAssay(mbrain) <- "RNA"
mbrain <- Seurat::SCTransform(mbrain)
mbrain <- Seurat::FindVariableFeatures(mbrain)
mbrain <- Seurat::RunPCA(mbrain, verbose = FALSE)
set.seed(10)
mbrain <- Seurat::RunUMAP(mbrain, dims = 1:50)

# --------------------------------------------------
# Get transferred cell labels from Jane.
# --------------------------------------------------

mbrain.jane <- readRDS("~/project/tiltedCCA/data/10x_mouseembryo/data_tenx_labels_Jane.rds")
mbrain$label_Savercat <- mbrain.jane$savercatLable

# --------------------------------------------------
# Cluster cells on basis of their scATAC-seq profiles
# --------------------------------------------------
set.seed(10)
DefaultAssay(mbrain) <- "ATAC"
mbrain <- Signac::RunTFIDF(mbrain)
mbrain <-  Signac::FindTopFeatures(mbrain, min.cutoff="q10")
mbrain <-  Signac::RunSVD(mbrain)  # question: is this svd run with only the top variable features?
set.seed(10)
mbrain <- Seurat::RunUMAP(mbrain, 
                          reduction="lsi", 
                          dims=2:50,
                          reduction.name="umap.atac", 
                          reduction.key="atacUMAP_")

# --------------------------------------------------
# Calculate a WNN graph
# --------------------------------------------------
set.seed(10)
mbrain <- Seurat::FindMultiModalNeighbors(mbrain, 
                                          reduction.list = list("pca", "lsi"), 
                                          dims.list = list(1:50, 2:50))
set.seed(10)
mbrain <- Seurat::RunUMAP(mbrain, 
                          nn.name = "weighted.nn", 
                          reduction.name = "umap.wnn", 
                          reduction.key = "wnnUMAP_")

Seurat::DefaultAssay(mbrain) <- "SCT"
mbrain[["percent.mt"]] <- Seurat::PercentageFeatureSet(object = mbrain, pattern = "^MT-")
mbrain[["percent.rb"]] <- Seurat::PercentageFeatureSet(object = mbrain, pattern = "^RPS")
mbrain <- Seurat::CellCycleScoring(mbrain, 
                                   g2m.features = cc.genes$g2m.genes, 
                                   s.features = cc.genes$s.genes)

save(mbrain, date_of_run, session_info,
     file = "../../../out/main/10x_mouseembryo_preprocessed.RData")

####################

plot1 <- Seurat::DimPlot(mbrain, reduction = "umap",
                         group.by = "label_Savercat", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Mouse Embryo E18 (10x, RNA+ATAC)\nRNA UMAP"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_mouseembryo_rna-umap.png"),
                plot1, device = "png", width = 8, height = 5, units = "in")

plot2 <- Seurat::DimPlot(mbrain, reduction = "umap.atac",
                         group.by = "label_Savercat", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot2 <- plot2 + ggplot2::ggtitle(paste0("Mouse Embryo E18 (10x, RNA+ATAC)\nATAC UMAP"))
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_mouseembryo_atac-umap.png"),
                plot2, device = "png", width = 8, height = 5, units = "in")

plot3 <- Seurat::DimPlot(mbrain, reduction = "umap.wnn",
                         group.by = "label_Savercat", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot3 <- plot3 + ggplot2::ggtitle(paste0("Mouse Embryo E18 (10x, RNA+ATAC)\nWNN UMAP"))
plot3 <- plot3 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_mouseembryo_wnn-umap.png"),
                plot3, device = "png", width = 8, height = 5, units = "in")


