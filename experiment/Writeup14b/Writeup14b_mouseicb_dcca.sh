#!/bin/bash
#$ -N mouseicb_dcca
#$ -j y
#$ -o ../../../../out/Writeup14b/qsub/
#$ -l m_mem_free=200G

Rscript --no-save Writeup14b_mouseicb_dcca.R
