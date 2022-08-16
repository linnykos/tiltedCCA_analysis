rm(list=ls())
library(Seurat)
library(MOFA2)

load("../../../out/main/citeseq_bm25_preprocessed.RData")
source("bm_25antibody_colorPalette.R")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

####

## https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/getting_started_R.html
## https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/downstream_analysis.html
MOFAobject <- MOFA2::create_mofa_from_Seurat(bm, 
                                             assays = c("RNA", "ADT"),
                                             slot = "scale.data")

data_opts <- MOFA2::get_default_data_options(MOFAobject)
head(data_opts)
model_opts <- MOFA2::get_default_model_options(MOFAobject)
head(model_opts)
train_opts <- MOFA2::get_default_training_options(MOFAobject)
head(train_opts)

MOFAobject <- MOFA2::prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

set.seed(10)
outfile <- "../../../out/main/citeseq_bm25_MOFA_model.hdf5"
MOFAobject.trained <- MOFA2::run_mofa(MOFAobject, outfile)

save(MOFAobject.trained, bm, 
     date_of_run, session_info,
     file = "../../../out/main/citeseq_bm25_MOFA.RData")