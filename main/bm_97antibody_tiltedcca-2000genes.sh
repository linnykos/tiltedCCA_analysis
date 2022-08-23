#!/bin/bash
#$ -N bm97_tcca-2000genes
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=60G

Rscript --no-save bm_97antibody_tiltedcca-2000genes.R
