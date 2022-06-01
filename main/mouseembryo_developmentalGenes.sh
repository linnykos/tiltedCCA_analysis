#!/bin/bash
#$ -N mouseembryo_developmentalGenes
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=100G

Rscript --no-save mouseembryo_developmentalGenes.R
