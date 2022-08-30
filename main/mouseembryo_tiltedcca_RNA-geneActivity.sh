#!/bin/bash
#$ -N mbrain_tcca-rna-gact
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=20G

Rscript --no-save mouseembryo_tiltedcca_RNA-geneActivity.R
