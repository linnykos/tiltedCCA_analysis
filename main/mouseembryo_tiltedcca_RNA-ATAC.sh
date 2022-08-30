#!/bin/bash
#$ -N mbrain_tcca-rna-atac
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=20G

Rscript --no-save mouseembryo_tiltedcca_RNA-ATAC.R
