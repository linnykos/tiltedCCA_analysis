#!/bin/bash
#$ -N mouseembryo_tcca
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=100G

Rscript --no-save mouseembryo_tiltedcca.R
