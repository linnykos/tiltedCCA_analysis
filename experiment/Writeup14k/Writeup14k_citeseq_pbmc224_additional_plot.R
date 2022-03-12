rm(list=ls())

library(Seurat)
library(SeuratData)

load("../../../../out/Writeup14h/Writeup14h_citeseq_pbmc224_preprocessed.RData")

cd_idx <- grep("^CD[0-9]", rownames(pbmc[["SCT"]]@scale.data))
rna_names <- sort(rownames(pbmc[["SCT"]]@scale.data)[cd_idx])
cd_idx <- grep("^CD[0-9]", rownames(pbmc[["ADT"]]@counts))
adt_names <- sort(rownames(pbmc[["ADT"]]@counts)[cd_idx])

idx <- intersect(rna_names, adt_names)
rna_adt_matching <- cbind(idx, idx)
colnames(rna_adt_matching) <- c("rna", "adt")
rna_adt_matching <- rbind(
  rna_adt_matching,
  matrix(
    c("CD4", "CD4-1",
      "CD44", "CD44-1",
      "CD38", "CD38-1"), byrow = T, ncol = 2
  ))

adt_names_nomatch <- adt_names[which(!adt_names %in% rna_names)]
adt_names_nomatch[which(adt_names_nomatch == "CD11a/CD18")] <- "CD11a"
adt_names_nomatch[which(adt_names_nomatch == "CD307c/FcRL3")] <- "CD307c"
adt_names_nomatch[which(adt_names_nomatch == "CD66a/c/e")] <- "CD66a"
adt_names_nomatch2 <- sapply(adt_names_nomatch, function(str){
  strsplit(str, split = "-")[[1]][1]
})
adt_names_nomatch2 <- unique(adt_names_nomatch2)

adt_synonyms <- Seurat::GeneSymbolThesarus(symbols = adt_names_nomatch2)
in_idx <- which(adt_synonyms %in% rownames(pbmc[["SCT"]]@scale.data))
rna_adt_matching <- rbind(
  rna_adt_matching,
  cbind(
    adt_synonyms[in_idx], names(adt_synonyms)[in_idx]
))

which(!rna_adt_matching[,"adt"] %in% rownames(pbmc[["ADT"]]@counts))
sort(rownames(pbmc[["ADT"]]@counts))

rna_adt_matching[c("CD11a","CD26","CD307c","CD56","CD66a"),"adt"] <- c("CD11a/CD18", "CD26-1", "CD307c/FcRL3", "CD56-1", "CD66a/c/e")
rownames(rna_adt_matching) <- NULL
rna_adt_matching <- as.data.frame(rna_adt_matching)
save(rna_adt_matching, file = "../../../../out/Writeup14h/Writeup14k_citeseq_pbmc224_rna-adt-matching.RData")

length(unique(rna_adt_matching[,1])) == nrow(rna_adt_matching)

##################################

# compute the correlation between RNA and ADT
corr_vec <- sapply(1:nrow(rna_adt_matching), function(i){
  rna_idx <- which(rownames(pbmc[["SCT"]]@scale.data) == rna_adt_matching[i,"rna"])
  adt_idx <- which(rownames(pbmc[["ADT"]]@scale.data) == rna_adt_matching[i,"adt"])
  
  stats::cor(pbmc[["SCT"]]@scale.data[rna_idx,], 
             pbmc[["ADT"]]@scale.data[adt_idx,])
})
names(corr_vec) <- rna_adt_matching[,"adt"]
quantile(corr_vec)
round(sort(corr_vec), 2)
round(corr_vec[order(names(corr_vec))], 2)

##################################

load("../../../../out/Writeup14j/Writeup14j_citeseq_pbmc224_dcca2.RData")
load("../../../../out/Writeup14i/Writeup14i_citeseq_pbmc224_dcca_variable_full_de.RData")
source("../Writeup14f/gene_exploration.R")

dcca_decomp <- multiomicCCA::dcca_decomposition(dcca_res)

summary_mat <- compute_variable_summary(mat = mat_2b, 
                                        common_mat = dcca_decomp$common_mat_2,
                                        metacell_clustering = NA)
colnames(summary_mat) <- c("r_squared", "p_value")
for(i in 1:nrow(summary_mat)){
  gene_name <- rownames(summary_mat)[i]
  
  max_celltype_idx <- max(combn_mat)
  celltype_pval <- sapply(1:max_celltype_idx, function(celltype){
    idx <- which(combn_mat == celltype, arr.ind = T)[,2]
    stats::quantile(sapply(idx, function(j){
      idx <- which(rownames(de_list[[j]]) == gene_name)
      if(length(idx) == 0) return(1)
      de_list[[j]][idx, "p_val"]
    }), probs = 0.75)
  })
  
  summary_mat[i,"p_value"] <- min(celltype_pval)
}
summary_mat[,"p_value"] <- -log10(summary_mat[,"p_value"])

summary_mat2 <- summary_mat[rownames(summary_mat) %in% rna_adt_matching[,"adt"],]
summary_cor <- sapply(rownames(summary_mat2), function(str){
  idx <- which(rna_adt_matching[,"adt"] == str)
  corr_vec[idx]
})
summary_mat2 <- cbind(summary_mat2, summary_cor)
colnames(summary_mat2)[3] <- "gene_adt_corr"
round(summary_mat2[order(summary_mat2[,"p_value"]),2:3],2)
round(summary_mat2[order(rownames(summary_mat2)),2:3],2)

x_idx <- which(rownames(pbmc[["SCT"]]@scale.data) == rna_adt_matching[which(rna_adt_matching[,"adt"] == "CD14"),"rna"])
y_idx <- which(rownames(pbmc[["SCT"]]@scale.data) == rna_adt_matching[which(rna_adt_matching[,"adt"] == "CD19"),"rna"])
x_vec <- pbmc[["SCT"]]@scale.data[x_idx,]
y_vec <- pbmc[["SCT"]]@scale.data[y_idx,]
label_idx <- as.numeric(pbmc@meta.data[,"celltype.l2"])
col_vec <- scales::hue_pal()(length(unique(label_idx)))
set.seed(10); shuf_idx <- sample(1:length(label_idx))
png("../../../../out/figures/Writeup14k/Writeup14k_biaxial_CD14_CD19_gene.png",
    height = 1500, width = 1500, res = 300, units = "px")
par(mar = c(4,4,0.5,0.5))
plot(x_vec[shuf_idx], y_vec[shuf_idx], 
     col = col_vec[label_idx[shuf_idx]],
     pch = 16,
     xlab = "CD14 (gene)", ylab = "CD19 (gene)")
graphics.off()
x_idx <- which(rownames(pbmc[["ADT"]]@scale.data) == rna_adt_matching[which(rna_adt_matching[,"adt"] == "CD14"),"adt"])
y_idx <- which(rownames(pbmc[["ADT"]]@scale.data) == rna_adt_matching[which(rna_adt_matching[,"adt"] == "CD19"),"adt"])
x_vec <- pbmc[["ADT"]]@scale.data[x_idx,]
y_vec <- pbmc[["ADT"]]@scale.data[y_idx,]
label_idx <- as.numeric(pbmc@meta.data[,"celltype.l2"])
col_vec <- scales::hue_pal()(length(unique(label_idx)))
set.seed(10); shuf_idx <- sample(1:length(label_idx))
png("../../../../out/figures/Writeup14k/Writeup14k_biaxial_CD14_CD19_protein.png",
    height = 1500, width = 1500, res = 300, units = "px")
par(mar = c(4,4,0.5,0.5))
plot(x_vec[shuf_idx], y_vec[shuf_idx], 
     col = col_vec[label_idx[shuf_idx]],
     pch = 16,
     xlab = "CD14 (antibody)", ylab = "CD19 (antibody)")
graphics.off()

##

x_idx <- which(rownames(pbmc[["SCT"]]@scale.data) == rna_adt_matching[which(rna_adt_matching[,"adt"] == "CD4-1"),"rna"])
y_idx <- which(rownames(pbmc[["SCT"]]@scale.data) == rna_adt_matching[which(rna_adt_matching[,"adt"] == "CD42b"),"rna"])
x_vec <- pbmc[["SCT"]]@scale.data[x_idx,]
y_vec <- pbmc[["SCT"]]@scale.data[y_idx,]
label_idx <- as.numeric(pbmc@meta.data[,"celltype.l2"])
col_vec <- scales::hue_pal()(length(unique(label_idx)))
set.seed(10); shuf_idx <- sample(1:length(label_idx))
png("../../../../out/figures/Writeup14k/Writeup14k_biaxial_CD4_CD42b_gene.png",
    height = 1500, width = 1500, res = 300, units = "px")
par(mar = c(4,4,0.5,0.5))
plot(x_vec[shuf_idx], y_vec[shuf_idx], 
     col = col_vec[label_idx[shuf_idx]],
     pch = 16,
     xlab = "CD4 (gene)", ylab = "GP1BA (gene)")
graphics.off()
x_idx <- which(rownames(pbmc[["ADT"]]@scale.data) == rna_adt_matching[which(rna_adt_matching[,"adt"] == "CD4-1"),"adt"])
y_idx <- which(rownames(pbmc[["ADT"]]@scale.data) == rna_adt_matching[which(rna_adt_matching[,"adt"] == "CD42b"),"adt"])
x_vec <- pbmc[["ADT"]]@scale.data[x_idx,]
y_vec <- pbmc[["ADT"]]@scale.data[y_idx,]
label_idx <- as.numeric(pbmc@meta.data[,"celltype.l2"])
col_vec <- scales::hue_pal()(length(unique(label_idx)))
set.seed(10); shuf_idx <- sample(1:length(label_idx))
png("../../../../out/figures/Writeup14k/Writeup14k_biaxial_CD4_CD42b_protein.png",
    height = 1500, width = 1500, res = 300, units = "px")
par(mar = c(4,4,0.5,0.5))
plot(x_vec[shuf_idx], y_vec[shuf_idx], 
     col = col_vec[label_idx[shuf_idx]],
     pch = 16,
     xlab = "CD4 (antibody)", ylab = "CD42b (antibody)")
graphics.off()

#################################

rm(list=ls())
library(Seurat)
load("../../../../out/Writeup14j/Writeup14j_citeseq_pbmc224_dcca2.RData")

print("Computing Consensus PCA")
n <- ncol(pbmc)
svd_1 <- svd(pbmc[["pca"]]@cell.embeddings)
embedding_1 <- multiomicCCA:::.mult_mat_vec(svd_1$u, svd_1$d/svd_1$d[1]*sqrt(n))
svd_2 <- svd(pbmc[["apca"]]@cell.embeddings)
embedding_2 <- multiomicCCA:::.mult_mat_vec(svd_2$u, svd_2$d/svd_2$d[1]*sqrt(n))
embedding_all <- cbind(embedding_1, embedding_2)
rownames(embedding_all) <- colnames(pbmc)
pca_res <- stats::prcomp(embedding_all, center = TRUE, scale. = FALSE)
consensus_mat <- pca_res$x[,1:min(rank_1, rank_2)] 

set.seed(10)
consensus_umap <- Seurat::RunUMAP(consensus_mat, 
                                  reduction.key = "umapConsensusPCA_")
pbmc[["consensus.umap"]] <- Seurat::CreateDimReducObject(consensus_umap@cell.embeddings,
                                                         assay = "ADT")

u_mat1 <- pbmc[["rna.umap"]]@cell.embeddings
u_mat2 <- pbmc[["consensus.umap"]]@cell.embeddings
tmp <- svd(t(u_mat1) %*% u_mat2)
rotation_mat <- tmp$u %*% t(tmp$v)
tmp <- u_mat2 %*% t(rotation_mat)
rownames(tmp) <- rownames(pbmc@meta.data)
colnames(tmp) <- colnames(pbmc[["consensus.umap"]]@cell.embeddings)
pbmc[["consensus.umap"]]@cell.embeddings <- tmp

plot1 <- Seurat::DimPlot(pbmc, reduction = "consensus.umap",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human PBMC (Cite-seq):\nConsensus PCA"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14k/Writeup14k_citeseq_pbmc224_dcca_consensuspca.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

Seurat::DefaultAssay(pbmc) <- "ADT"
plot1 <- Seurat::FeaturePlot(pbmc, 
                             features = "CD14",
                             slot = "scale.data",
                             cols = c("lightgrey", "darkgreen"),
                             reduction = "consensus.umap")
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human PBMC (Cite-seq):\nConsensus PCA, CD14 (antibody)"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14k/Writeup14k_citeseq_pbmc224_consensuspca_cd14_antibody.png"),
                plot1, device = "png", width = 4.5, height = 5, units = "in")
Seurat::DefaultAssay(pbmc) <- "SCT"
plot1 <- Seurat::FeaturePlot(pbmc, 
                             features = "CD14",
                             slot = "scale.data",
                             cols = c("lightgrey", "deepskyblue3"),
                             reduction = "consensus.umap")
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human PBMC (Cite-seq):\nConsensus PCA, CD14 (gene)"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14k/Writeup14k_citeseq_pbmc224_consensuspca_cd14_gene.png"),
                plot1, device = "png", width = 4.5, height = 5, units = "in")


Seurat::DefaultAssay(pbmc) <- "ADT"
plot1 <- Seurat::FeaturePlot(pbmc, 
                             features = "CD4-1",
                             slot = "scale.data",
                             cols = c("lightgrey", "darkgreen"),
                             reduction = "consensus.umap")
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human PBMC (Cite-seq):\nConsensus PCA, CD4 (antibody)"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14k/Writeup14k_citeseq_pbmc224_consensuspca_cd4_antibody.png"),
                plot1, device = "png", width = 4.5, height = 5, units = "in")
Seurat::DefaultAssay(pbmc) <- "SCT"
plot1 <- Seurat::FeaturePlot(pbmc, 
                             features = "CD4",
                             slot = "scale.data",
                             cols = c("lightgrey", "deepskyblue3"),
                             reduction = "consensus.umap")
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human PBMC (Cite-seq):\nConsensus PCA, CD4 (gene)"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14k/Writeup14k_citeseq_pbmc224_consensuspca_cd4_gene.png"),
                plot1, device = "png", width = 4.5, height = 5, units = "in")


Seurat::DefaultAssay(pbmc) <- "ADT"
plot1 <- Seurat::FeaturePlot(pbmc, 
                             features = "CD69",
                             slot = "scale.data",
                             cols = c("lightgrey", "darkgreen"),
                             reduction = "consensus.umap")
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human PBMC (Cite-seq):\nConsensus PCA, CD69 (antibody)"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14k/Writeup14k_citeseq_pbmc224_consensuspca_cd69_antibody.png"),
                plot1, device = "png", width = 4.5, height = 5, units = "in")
Seurat::DefaultAssay(pbmc) <- "SCT"
plot1 <- Seurat::FeaturePlot(pbmc, 
                             features = "CD69",
                             slot = "scale.data",
                             cols = c("lightgrey", "deepskyblue3"),
                             reduction = "consensus.umap")
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human PBMC (Cite-seq):\nConsensus PCA, CD69 (gene)"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14k/Writeup14k_citeseq_pbmc224_consensuspca_cd69_gene.png"),
                plot1, device = "png", width = 4.5, height = 5, units = "in")



##################################

# for fun, let's look at the 25 antibodies
bm <- SeuratData::LoadData(ds = "bmcite")
sort(rownames(bm[["ADT"]]))