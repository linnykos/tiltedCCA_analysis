#!/bin/bash
#$ -N abseq_DE-distinct
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=100G

Rscript --no-save abseq_bm_97antibodyRef_differential-distinct.R
