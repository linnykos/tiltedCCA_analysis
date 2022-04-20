#!/bin/bash
#$ -N bm_25antibody_differential
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=50G

Rscript --no-save bm_25antibody_differential.R
