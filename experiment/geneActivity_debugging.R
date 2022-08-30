zz <- mbrain[["ATAC"]]
zz@fragments[[1]]
head(zz@fragments[[1]]@cells)
tmp <- zz@fragments[[1]]@cells
cellnames <- zz@fragments[[1]]@cells
names(cellnames) <- NULL

head(colnames(mbrain[["RNA"]]@counts))
all(tmp %in% colnames(mbrain[["RNA"]]@counts))

# https://github.com/timoast/signac/issues/205
mbrain@assays[["ATAC"]]@fragments[[1]]@cells %>% head
colnames(x = mbrain[["ATAC"]]) %>% head

assay <- "ATAC"
counts <- FeatureMatrix(
  fragments = mbrain@assays[["ATAC"]]@fragments[1],
  features = transcripts,
  cells = cellnames) ## return all zero matrix

frag.file <- file.path("~/nzhanglab/data/10x_mouse_embryo/e18_mouse_brain_fresh_5k_atac_fragments.tsv.gz")
fragment_file <- Signac::CreateFragmentObject(frag.file)

############

# annotation <- Annotation(mbrain)
annotation <- annotations
transcripts <- Signac:::CollapseToLongestTranscript(ranges = annotation)
transcripts <- transcripts[transcripts$gene_biotype == "protein_coding"]
transcripts <- Signac:::Extend(x = transcripts, upstream = 500, 
                      downstream = 500)
counts <- Signac::FeatureMatrix(
  fragments = fragment_file,
  features = transcripts,
  cells = cellnames)

#####################3

fragment <- mbrain@assays[["ATAC"]]@fragments[[1]]
frag.cells <- Signac::GetFragmentData(object = fragment, slot = "cells")
