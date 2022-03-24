#!/bin/bash
#$ -N greenleaf_atac
#$ -j y
#$ -o ../../../../out/Writeup14m/qsub/
#$ -l m_mem_free=200G

Rscript --no-save Writeup14m_10x_greenleaf_ataconly_preprocess.R
