#!/bin/bash
#$ -N mouseembryo_differential
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=50G

Rscript --no-save mouseembryo_differential.R
