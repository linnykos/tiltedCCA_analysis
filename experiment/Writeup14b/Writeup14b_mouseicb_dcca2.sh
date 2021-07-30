#!/bin/bash
#$ -N mouseicb_dcca2
#$ -j y
#$ -o ../../../../out/Writeup14b/qsub/
#$ -l m_mem_free=150G

Rscript --no-save Writeup14b_mouseicb_dcca2.R
