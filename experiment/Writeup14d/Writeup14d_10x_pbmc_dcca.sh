#!/bin/bash
#$ -N 10x_pbmc_dcca
#$ -j y
#$ -o ../../../../out/Writeup14d/qsub/
#$ -l m_mem_free=75G

Rscript --no-save Writeup14d_10x_pbmc_dcca.R
