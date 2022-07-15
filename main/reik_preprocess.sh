#!/bin/bash
#$ -N reik_preprocess
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=250G

Rscript --no-save reik_preprocess.R
