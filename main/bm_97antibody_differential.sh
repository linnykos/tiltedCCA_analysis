#!/bin/bash
#$ -N bm97_differential
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=100G

Rscript --no-save bm_97antibody_differential.R
