#!/bin/bash
#$ -N pbmc224_tcca-25ab
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=100G

Rscript --no-save pbmc_224antibody_tiltedcca-25antibodies.R
