#!/bin/bash
#$ -N mouseembryo_dcca
#$ -j y
#$ -o ../../../../out/Writeup14f/qsub/
#$ -l m_mem_free=50G

Rscript --no-save Writeup14f_10x_mouseembryo_dcca.R
