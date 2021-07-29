rm(list=ls())
load("../../../../out/Writeup14b/Writeup14b_citeseq_pbmc228_dcca1.RData")
tmp <- dcca_res$common_score; colnames(tmp) <- paste0("p", 1:ncol(tmp))
pbmc2 <- Seurat::CreateSeuratObject(counts = t(tmp))
pbmc2[["celltype"]] <- membership_vec
rm(list = c("dcca_res", "tmp")); gc(T)

load("../../../../out/Writeup14b/Writeup14b_citeseq_pbmc228_dcca2.RData")
pbmc2[["common"]] <- Seurat::CreateDimReducObject(embedding = rna_embeddings[[1]], key = "UMAP", assay = "RNA")
pbmc2[["distinct"]] <- Seurat::CreateDimReducObject(embedding = rna_embeddings[[2]], key = "UMAP", assay = "RNA")
pbmc2[["everything"]] <- Seurat::CreateDimReducObject(embedding = rna_embeddings[[3]], key = "UMAP", assay = "RNA")
pbmc2[["clap"]] <- Seurat::CreateDimReducObject(embedding = c_eig, key = "clap", assay = "RNA")
pbmc2[["dlap"]] <- Seurat::CreateDimReducObject(embedding = d_eig, key = "dlap", assay = "RNA")
pbmc2[["elap"]] <- Seurat::CreateDimReducObject(embedding = e_eig, key = "elap", assay = "RNA")
rm(list = c("c_eig", "d_eig", "e_eig", "rna_embeddings", "rna_frnn")); gc(T)

load("../../../../out/Writeup14b/Writeup14b_citeseq_pbmc228_dcca3.RData")
pbmc2[["common2"]] <- Seurat::CreateDimReducObject(embedding = protein_embeddings[[1]], key = "UMAP", assay = "RNA")
pbmc2[["distinct2"]] <- Seurat::CreateDimReducObject(embedding = protein_embeddings[[2]], key = "UMAP", assay = "RNA")
pbmc2[["everything2"]] <- Seurat::CreateDimReducObject(embedding = protein_embeddings[[3]], key = "UMAP", assay = "RNA")
pbmc2[["clap2"]] <- Seurat::CreateDimReducObject(embedding = c_eig2, key = "clap", assay = "RNA")
pbmc2[["dlap2"]] <- Seurat::CreateDimReducObject(embedding = d_eig2, key = "dlap", assay = "RNA")
pbmc2[["elap2"]] <- Seurat::CreateDimReducObject(embedding = e_eig2, key = "elap", assay = "RNA")
rm(list = c("c_eig2", "d_eig2", "e_eig2", "protein_embeddings", "protein_frnn")); gc(T)

load("../../../../out/Writeup14b/Writeup14b_citeseq_pbmc228_dcca4.RData")
pbmc2[["combined"]] <- Seurat::CreateDimReducObject(embedding = combined_common_umap, key = "UMAP", assay = "RNA")
rm(list = c("combined_common_umap")); gc(T)

load("../../../../out/Writeup14b/Writeup14b_citeseq_pbmc228_dcca5.RData")
pbmc2[["both"]] <- Seurat::CreateDimReducObject(embedding = both_embeddings$everything, key = "UMAP", assay = "RNA")
rm(list = c("both_embeddings")); gc(T)

###########################

# done with all the calculations (for now). now plot

title_vec <- c("Common view", "Distinct view", "Everything view")
main_vec <- c("common", "distinct", "everything")

# plot RNA embeddings
for(i in 1:3){
  plot1 <- Seurat::DimPlot(pbmc2, reduction = main_vec[i],
                           group.by = "celltype", label = TRUE,
                           repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("CITE-Seq PBMC-228 (RNA)\n", title_vec[i]))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_citeseq_pbmc228_dcca_rna_", main_vec[i], "_umap.png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
}

# plot ATAC embeddings
for(i in 1:3){
  plot1 <- Seurat::DimPlot(pbmc2, reduction = paste0(main_vec[i], "2"),
                           group.by = "celltype", label = TRUE,
                           repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("CITE-Seq PBMC-228 (Protein)\n", title_vec[i]))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_citeseq_pbmc228_dcca_atac_", main_vec[i], "_umap.png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
}

# plot the combined view
plot1 <- Seurat::DimPlot(pbmc2, reduction = "combined",
                         group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle("CITE-Seq PBMC-228\nBoth Common")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_citeseq_pbmc228_dcca_both_common_umap.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

# plot the everything combined for both
plot1 <- Seurat::DimPlot(pbmc2, reduction = "both",
                         group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle("CITE-Seq PBMC-228\nBoth Everything")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_citeseq_pbmc228_dcca_both_everything_umap.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")


#######################################
# now plot the bases -- first RNA
plot1 <- Seurat::FeaturePlot(pbmc2, features = paste0("clap_", 1:16), reduction = "combined")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_citeseq_pbmc228_dcca_rna_common_basis.png"),
                plot1, device = "png", width = 16, height = 12, units = "in")

plot1 <- Seurat::FeaturePlot(pbmc2, features = paste0("dlap_", 1:16), reduction = "distinct")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_citeseq_pbmc228_dcca_rna_distinct_basis.png"),
                plot1, device = "png", width = 16, height = 12, units = "in")

plot1 <- Seurat::FeaturePlot(pbmc2, features = paste0("elap_", 1:16), reduction = "both")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_citeseq_pbmc228_dcca_rna_everything_basis.png"),
                plot1, device = "png", width = 16, height = 12, units = "in")

# next ATAC
plot1 <- Seurat::FeaturePlot(pbmc2, features = paste0("clap2_", 1:16), reduction = "combined")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_citeseq_pbmc228_dcca_atac_common_basis.png"),
                plot1, device = "png", width = 16, height = 12, units = "in")

plot1 <- Seurat::FeaturePlot(pbmc2, features = paste0("dlap2_", 1:16), reduction = "distinct2")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_citeseq_pbmc228_dcca_atac_distinct_basis.png"),
                plot1, device = "png", width = 16, height = 12, units = "in")

plot1 <- Seurat::FeaturePlot(pbmc2, features = paste0("elap2_", 1:16), reduction = "both")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_citeseq_pbmc228_dcca_atac_everything_basis.png"),
                plot1, device = "png", width = 16, height = 12, units = "in")

##################
rm(list = c("pbmc2")); gc(T)

nn <- 30

# compute local enrichment
load("../../../../out/Writeup14b/Writeup14b_citeseq_pbmc228_dcca2.RData")
rm(list = c("c_eig", "d_eig", "e_eig", "rna_embeddings")); gc(T)
set.seed(10)
rna_local <- multiomicCCA::clisi_information(rna_frnn$c_g, rna_frnn$d_g, rna_frnn$e_g, 
                                             membership_vec = membership_vec)
rm(list = c("rna_frnn")); gc(T)
# rna_local$common_clisi$membership_info
# rna_local$distinct_clisi$membership_info

load("../../../../out/Writeup14b/Writeup14b_citeseq_pbmc228_dcca3.RData")
rm(list = c("c_eig2", "d_eig2", "e_eig2", "protein_embeddings")); gc(T)
set.seed(10)
protein_local <- multiomicCCA::clisi_information(protein_frnn$c_g, protein_frnn$d_g, protein_frnn$e_g, 
                                                 membership_vec = membership_vec)
rm(list = c("protein_frnn")); gc(T)
# atac_local$common_clisi$membership_info
# atac_local$distinct_clisi$membership_info

tmp <- multiomicCCA::plot_clisi(rna_local, protein_local)
tmp2 <- cowplot::plot_grid(tmp[[1]], tmp[[2]])
cowplot::save_plot(filename = "../../../../out/figures/Writeup14b/Writeup14b_citeseq_pbmc228_enrichment.png", 
                   tmp2, ncol = 1, nrow = 2, base_height = 1.75, base_asp = 4, device = "png")

