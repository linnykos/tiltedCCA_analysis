#!/bin/bash
#$ -N reik_differential
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=150G

Rscript --no-save reik_differential.R
