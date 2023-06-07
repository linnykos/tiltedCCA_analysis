#!/bin/bash
#$ -N bm_tcca_noise-0.25
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=100G

Rscript --no-save bm_25antibody_tiltedcca_added-noise_0.25.R
