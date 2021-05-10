rm(list=ls())
library(simulator)
library(multiomicCCA)
source("../multiomicCCA_analysis/simulation/data_generator.R")

df_param <- data.frame(setting = c(1, 2, 3, 4, 5, 6), 
                       rank_1 =  c(2, 2, 2, 2, 2, 2),
                       rank_2 =  c(2, 2, 2, 2, 2, 2))
vec <- df_param[5,]
set.seed(10)
dat <- simulation_all(vec$setting)

mofa <- MOFA2::create_mofa_from_matrix(list(Mode1 = t(dat$dat$mat_1), Mode2 = t(dat$dat$mat_2)))
model_opts <- MOFA2::get_default_model_options(mofa)
model_opts$num_factors <- 2
mofa <- MOFA2::prepare_mofa(mofa, model_options = model_opts)
mofa <- MOFA2::run_mofa(mofa)

MOFA2::plot_factor_cor(mofa)
MOFA2::plot_variance_explained(mofa)
MOFA2::plot_variance_explained(mofa, plot_total = TRUE)[[2]]
# correlate_factors_with_covariates(mofa, covariates = c("nFeature_RNA","nFeature_ATAC"))
mofa@samples_metadata[,"celltype"] <- dat$true_membership_vec
MOFA2::plot_factor(mofa, factors=1, group_by = "celltype", color_by="celltype") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(color="black", angle=40, vjust=1, hjust=1))
MOFA2::plot_factor(mofa, factors=2, group_by = "celltype", color_by="celltype") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(color="black", angle=40, vjust=1, hjust=1))
MOFA2::plot_dimred(mofa, method = "UMAP", color_by = "celltype", label = TRUE, 
            stroke=0.05, dot_size = 1, legend = FALSE) 

###################

set.seed(10)
scAI_outs <- scAI::create_scAIobject(raw.data = list(Mode1 = t(dat$dat$mat_1), Mode2 = t(dat$dat$mat_2)), do.sparse = F)
scAI_outs <- scAI::preprocessing(scAI_outs, assay = NULL, minFeatures = 0, minCells = 0, libararyflag = F, logNormalize = F)
scAI_outs <- scAI::run_scAI(scAI_outs, K = 5, nrun = 5, s = 1)
scAI_outs <- scAI::addpData(scAI_outs, pdata = as.data.frame(dat$true_membership_vec), pdata.name = "labels")
scAI::lmHeatmap(scAI_outs, color.by = "labels")

scAI_outs <- scAI::reducedDims(scAI_outs, method = "umap")
scAI::cellVisualization(scAI_outs, scAI_outs@embed$umap, color.by = "labels")
