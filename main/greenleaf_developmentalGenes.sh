#!/bin/bash
#$ -N greenleaf_developmentalGenes
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=100G

Rscript --no-save greenleaf_developmentalGenes.R
