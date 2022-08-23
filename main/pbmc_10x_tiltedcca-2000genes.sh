#!/bin/bash
#$ -N pbmc_10x_tcca-2000genes
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=100G

Rscript --no-save pbmc_10x_tiltedcca-2000genes.R
