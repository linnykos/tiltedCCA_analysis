# from https://satijalab.org/signac/articles/mouse_brain_vignette.html

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(patchwork)
set.seed(1234)

counts <- Read10X_h5("~/nzhanglab/data/10x_mouse_embryo_onlyATAC/atac_v1_adult_brain_fresh_5k_filtered_peak_bc_matrix.h5")
metadata <- read.csv(
  file = "~/nzhanglab/data/10x_mouse_embryo_onlyATAC/atac_v1_adult_brain_fresh_5k_singlecell.csv",
  header = TRUE,
  row.names = 1
)

brain_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = "mm10",
  fragments = '~/nzhanglab/data/10x_mouse_embryo_onlyATAC/atac_v1_adult_brain_fresh_5k_fragments.tsv.gz',
  min.cells = 1
)

brain <- CreateSeuratObject(
  counts = brain_assay,
  assay = 'peaks',
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

############################

fragment <- brain@assays[["peaks"]]@fragments[[1]]
frag.cells <- Signac::GetFragmentData(object = fragment, slot = "cells")

