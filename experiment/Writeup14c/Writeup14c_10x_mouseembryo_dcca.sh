#!/bin/bash
#$ -N dcca_mouse
#$ -j y
#$ -o ../../../../out/Writeup14c/qsub/
#$ -l m_mem_free=50G

Rscript --no-save Writeup14c_10x_mouseembryo_dcca.R
