#!/bin/bash
#$ -N reik_tiltedcca2
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=100G

Rscript --no-save reik_tiltedcca2.R
