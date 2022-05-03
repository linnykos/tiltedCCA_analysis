#!/bin/bash
#$ -N pbmc_224antibody_sel_alt
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=100G

Rscript --no-save pbmc_224antibody_antibody_selection_alternatives.R
