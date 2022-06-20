#!/bin/bash
#$ -N pbmc224_antibody_distinct
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=100G

Rscript --no-save pbmc_224antibody_differential_distinct.R
