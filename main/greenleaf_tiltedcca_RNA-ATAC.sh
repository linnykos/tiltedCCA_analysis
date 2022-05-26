#!/bin/bash
#$ -N greenleaf_tcca
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=100G

Rscript --no-save greenleaf_tiltedcca_RNA-ATAC.R
