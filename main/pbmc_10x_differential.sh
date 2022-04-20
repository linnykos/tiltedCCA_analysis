#!/bin/bash
#$ -N pbmc_10x_differential
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=50G

Rscript --no-save pbmc_10x_differential.R
