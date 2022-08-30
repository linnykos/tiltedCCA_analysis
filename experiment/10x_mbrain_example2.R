# from https://satijalab.org/signac/articles/mouse_brain_vignette.html

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(patchwork)
set.seed(1234)

counts <- Read10X_h5("e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix.h5")
metadata <- read.csv(
  file = "e18_mouse_brain_fresh_5k_per_barcode_metrics.csv",
  header = TRUE,
  row.names = 1
)

brain_assay <- CreateChromatinAssay(
  counts =  counts[["Peaks"]],
  sep = c(":", "-"),
  genome = "mm10",
  fragments = 'e18_mouse_brain_fresh_5k_atac_fragments.tsv.gz',
  min.cells = 1
)

brain <- CreateSeuratObject(
  counts = brain_assay,
  assay = 'Peaks',
  project = 'ATAC',
  meta.data = metadata
)

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
# add the gene information to the object
Annotation(brain) <- annotations

gene.activities <- GeneActivity(brain)
length(gene.activities@x)

##########################

fragment <- brain[["Peaks"]]@fragments[[1]]
frag.cells <- Signac::GetFragmentData(object = fragment, slot = "cells")
head(frag.cells)
all(frag.cells == colnames(brain[["Peaks"]]))

Annotation(brain)
table(Annotation(brain)@seqnames)
table(Annotation(brain)@elementMetadata$type)

granges(brain)
table(granges(brain)@seqnames)

############################

fragment <- brain@assays[["peaks"]]@fragments[[1]]
frag.cells <- Signac::GetFragmentData(object = fragment, slot = "cells")

