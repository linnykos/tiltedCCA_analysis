rm(list=ls())

library(Seurat)
set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

load("../../../../data/TetTCR-SeqHD_Antigen_t_cells/T1D.RData")
ls()

DBEC[1:5,1:5]
dim(DBEC)
protein_idx <- grep("pAbO", colnames(DBEC))
rna_idx <- c(1:ncol(DBEC))[-protein_idx]
length(protein_idx); length(rna_idx)

T1D_meta[1:5,1:5]
dim(T1D_meta)
head(T1D_meta$Top1peptide)
length(unique(T1D_meta$Top1peptide))
table(T1D_meta$Top1peptide)
quantile(table(T1D_meta$Top1peptide))

head(T1D_meta$TCR)
length(unique(T1D_meta$TCR))
quantile(table(T1D_meta$TCR))

head(T1D_meta$TCRa_1st_cdr_aa)
length(unique(T1D_meta$TCRa_1st_cdr_aa))
quantile(table(T1D_meta$TCRa_1st_cdr_aa))

head(T1D_meta$TCRa_1st_cdr_aa)
head(T1D_meta$TCRb_1st_cdr_aa)
head(T1D_meta$TCR)

length(which(table(T1D_meta$TCR) > 10))
length(which(table(T1D_meta$TCR) > 10))/length(unique(T1D_meta$TCR))

table(T1D_meta$disease)
table(T1D_meta$donor)
table(T1D_meta$leiden)

#####################################

apply(T1D_meta, 2, function(x){length(unique(x))})

# let's make some plots
zz <- t(DBEC[,rna_idx])
metadata <- T1D_meta
metadata <- metadata[,which(colnames(metadata) %in% c("chip", "donor", "gene",  "Top1peptide", "TCR", "disease", "group", 
                                                      "leiden", "tcr.frequency", "gene.final", "thred15", "seurat_clusters",
                                                      "Ab_N", "combocluster"))]
tettcr <- Seurat::CreateSeuratObject(counts = zz, meta.data = metadata)
tettcr <- Seurat::SCTransform(tettcr)
length(tettcr[["SCT"]]@var.features)
tettcr <- Seurat::RunPCA(tettcr, verbose = F)
set.seed(10)
tettcr <- Seurat::RunUMAP(tettcr, dims = 1:50,
                          reduction.name = "umap.rna", reduction.key = "rnaUMAP_")

zz <- t(DBEC[,protein_idx])
tettcr[["ADT"]] <- Seurat::CreateAssayObject(counts = zz)
Seurat::DefaultAssay(tettcr) <- "ADT"
tettcr <- Seurat::NormalizeData(tettcr, normalization.method = 'CLR', margin = 2)
tettcr <- Seurat::ScaleData(tettcr)
Seurat::VariableFeatures(tettcr) <- rownames(tettcr[["ADT"]])
tettcr <- Seurat::RunPCA(tettcr, reduction.name = 'apca', verbose = F)
set.seed(10)
tettcr <- Seurat::RunUMAP(tettcr, dims = 1:30,
                          reduction.name = "umap.adt", reduction.key = "adtUMAP_")

set.seed(10)
tettcr <- Seurat::FindMultiModalNeighbors(
  tettcr, reduction.list = list("pca", "apca"), 
  dims.list = list(1:50, 1:30), modality.weight.name = "RNA.weight"
)
set.seed(10)
tettcr <- Seurat::RunUMAP(tettcr, nn.name = "weighted.nn", 
                          reduction.name = "umap.wnn", reduction.key = "wnnUMAP_")

zz <- as.matrix(T1D_meta[,c("UMAP_1", "UMAP_2")])
all(rownames(zz) == colnames(tettcr))
tettcr[["umap.totalvi"]] <- Seurat::CreateDimReducObject(embeddings = zz, assay = "RNA")

celltype_names <- c(paste0("TNaive_", 1:4), 
                    "NK-like",
                    "TCenMemory", "TTrans",
                    paste0("TEffMemory_", 1:3),
                    paste0("TEff", 1:2), "HCVSpikeIn")
celltype_vec <- as.character(as.numeric(tettcr$leiden))
for(i in 1:length(celltype_names)){
  celltype_vec[which(celltype_vec == i)] <- celltype_names[i]
}
tettcr$celltype <- celltype_vec

#######################

color_vec <- c(rgb(167, 205, 228, max = 255), #1
               rgb(43, 128, 183, max = 255), #2
               rgb(154, 205, 146, max = 255), #3
               rgb(83, 174, 68, max = 255), #4
               rgb(184, 154, 117, max = 255), #5
               rgb(237, 79, 80, max = 255), #6
               rgb(238, 109, 68, max = 255), #7
               rgb(253, 162, 63, max = 255), #8
               rgb(238, 144, 72, max = 255), #9
               rgb(179, 149, 201, max = 255), #10
               rgb(131, 93, 156, max = 255), #11
               rgb(249, 241, 142, max = 255), #12
               rgb(176, 90, 43, max = 255))
names(color_vec) <- celltype_names

save(tettcr, date_of_run, session_info, color_vec,
     file = "../../../../out/Writeup14/Writeup14_TetCRSeqHD_preprocessed.RData")

###########################

# make some plots
plot1 <- Seurat::DimPlot(tettcr, reduction = "umap.rna",
                         group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         cols = color_vec)
plot1 <- plot1 + ggplot2::ggtitle(paste0("TetCR-SeqHD\n (CD8+ T-cells, RNA)"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14h/TetCRSeqHD_rna_umap.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <- Seurat::DimPlot(tettcr, reduction = "umap.adt",
                         group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         cols = color_vec)
plot1 <- plot1 + ggplot2::ggtitle(paste0("TetCR-SeqHD\n (CD8+ T-cells, ADT)"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14h/TetCRSeqHD_adt_umap.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <- Seurat::DimPlot(tettcr, reduction = "umap.wnn",
                         group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         cols = color_vec)
plot1 <- plot1 + ggplot2::ggtitle(paste0("TetCR-SeqHD\n (CD8+ T-cells, WNN)"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14h/TetCRSeqHD_wnn_umap.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <- Seurat::DimPlot(tettcr, reduction = "umap.totalvi",
                         group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         cols = color_vec)
plot1 <- plot1 + ggplot2::ggtitle(paste0("TetCR-SeqHD\n (CD8+ T-cells, TotalVI)"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14h/TetCRSeqHD_totalvi_umap.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

##########################

# make a heatmap
rna_mat <- tettcr[["SCT"]]@scale.data
adt_mat <- tettcr[["ADT"]]@scale.data
protein_name <- c("CD28", "CD38", "CD56", "CD366", "CD26",
                  "CD45RA", "CD197", "CD62L", "CD25", "CD20",
                  "CD95", "CD45RO", "CD279", "CD57")
gene_name <- c("KLRB1", "GNLY", "ZNF683", "LGALS3", "GZMK",
               "DUSP4", "HLA-DPA1", "HLA-DRB3", "HLA-DRA",
               "GZMB", "ADGRG1", "FCGR3A", "FGFBP2", "GZMH",
               "ITGAM", "KLRG1", "IFNG", "CD244", "CD160",
               "ZEB2", "NKG7", "CCL4", "PRF1", "KLRD1",
               "S1PR5", "KLRF1")
# subset relevant variables
all_mat <- rbind(rna_mat[sapply(paste0(gene_name, "-"), function(x){grep(x, rownames(rna_mat))[1]}),],
                 adt_mat[sapply(protein_name, function(x){grep(x, rownames(adt_mat))}),])
# convert to normalized z-scores
for(j in 1:nrow(all_mat)){
  mean_val <- mean(all_mat[j,])
  sd_val <- sd(all_mat[j,])
  all_mat[j,] <- c(all_mat[j,] - mean_val)/sd_val
}
all_mat <- t(scale(t(all_mat), center = T, scale = T))
# average over celltypes
all_mat2 <- sapply(1:12, function(x){
  idx <- which(tettcr$leiden == x)
  matrixStats::rowMeans2(all_mat[,idx])
})
rownames(all_mat2) <- c(paste0("gene-", gene_name), paste0("protein-", protein_name))
colnames(all_mat2) <- 1:12

png("../../../../out/figures/Writeup14h/TetCRSeqHD_heatmap.png",
    height = 2500, width = 1750, units = "px", res = 300)
gplots::heatmap.2(all_mat2, scale = "none", col = gplots::bluered(100),
                  trace = "none", density.info = "none")
graphics.off()

############################

grep("CAYRSPPSSEKLVF_CASSFLGTGLNEQYF", tettcr$TCR)
tcr_tab <- table(tettcr$TCR)
tcr_tab <- tcr_tab[-which(names(tcr_tab) == "_")]
idx <- which(sapply(names(tcr_tab), function(x){
  res <- strsplit(x, split = "_")[[1]]
  all(sapply(res, nchar) > 0) & length(res) == 2
}))
tcr_tab <- tcr_tab[idx]
idx <- order(tcr_tab, decreasing = T)[1:35]
tcr_names <- names(tcr_tab)[idx]
sum(tcr_tab[tcr_names])

tcr_idx <- which(tettcr$TCR %in% tcr_names)
antibody_tab <- table(tettcr$Top1peptide[tcr_idx])
idx <- order(antibody_tab, decreasing = T)[1:20]
antibody_names <- names(antibody_tab)[idx]

idx <- intersect(tcr_idx,
                 which(tettcr$Top1peptide %in% antibody_names))

tab_mat <- table(tettcr$Top1peptide[idx], tettcr$TCR[idx])
tab_mat <- log10(tab_mat+1)
tab_mat <- tab_mat[,order(matrixStats::colSums2(tab_mat), decreasing = T)]
row_idx <- unique(unlist(sapply(1:ncol(tab_mat), function(j){which(tab_mat[,j] > 0)})))
tab_mat <- tab_mat[row_idx,]

.rotate <- function(mat){t(mat)[,nrow(mat):1]}

png("../../../../out/figures/Writeup14h/TetCRSeqHD_tcr_antigen.png",
    height = 2250, width = 3500, units = "px", res = 300)
par(mar = c(15,0.5,0.5,10))
image(.rotate(tab_mat), col = grDevices::colorRampPalette(c("white", "red"))(10),
      xaxt = "n", yaxt = "n")

tmp <- 1/((ncol(tab_mat)-1)*2)
x_axis <- seq(0, 1, by = 2*tmp)-tmp
text(x = x_axis, y = par("usr")[3]-0.05,
     labels = colnames(tab_mat),
     xpd = NA, pos = 4, srt = -90, cex = 0.75)

tmp <- 1/((nrow(tab_mat)-1)*2)
y_axis <- seq(1, 0, by = -2*tmp)
text(x = par("usr")[4], y = y_axis,
     labels = rownames(tab_mat),
     xpd = NA, pos = 4, cex = 0.75)
graphics.off()

