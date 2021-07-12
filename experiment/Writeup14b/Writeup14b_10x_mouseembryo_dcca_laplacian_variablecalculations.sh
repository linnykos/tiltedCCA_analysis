#!/bin/bash
#$ -N dcca_mouseembryo
#$ -j y
#$ -o ../../../../out/Writeup14b/qsub/
#$ -l m_mem_free=30G

Rscript --no-save Writeup14b_10x_mouseembryo_dcca_laplacian_variablecalculations.R
