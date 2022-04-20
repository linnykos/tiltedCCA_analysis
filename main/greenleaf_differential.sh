#!/bin/bash
#$ -N greenleaf_differential
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=50G

Rscript --no-save greenleaf_differential.R
