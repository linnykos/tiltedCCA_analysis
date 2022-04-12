#!/bin/bash
#$ -N bm_tcca
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=100G

Rscript --no-save bm_25antibody_tiltedcca_simplified.R
