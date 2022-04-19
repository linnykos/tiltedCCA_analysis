#!/bin/bash
#$ -N pbmc_10x_tcca
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=100G

Rscript --no-save pbmc_10x_tiltedcca.R
