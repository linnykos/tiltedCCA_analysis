#!/bin/bash
#$ -N pbmc_224antibody_tcca
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=100G

Rscript --no-save pbmc_224antibody_tiltedcca.R
