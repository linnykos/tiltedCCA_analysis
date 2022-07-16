#!/bin/bash
#$ -N reik_tiltedcca
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=200G

Rscript --no-save reik_tiltedcca.R
