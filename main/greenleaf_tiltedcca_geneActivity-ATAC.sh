#!/bin/bash
#$ -N greenleaf_tcca_GA-ATAC
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=100G

Rscript --no-save greenleaf_tiltedcca_geneActivity-ATAC.R
