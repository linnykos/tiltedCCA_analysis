#!/bin/bash
#$ -N bm_25antibody_scai
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=100G

Rscript --no-save bm_25antibody_scai.R
