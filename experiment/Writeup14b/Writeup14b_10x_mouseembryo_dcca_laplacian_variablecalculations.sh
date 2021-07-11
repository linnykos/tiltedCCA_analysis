#!/bin/bash
#$ -N bm_glmpca_nb
#$ -j y
#$ -o ../../../../out/Writeup14b/qsub/

Rscript --no-save Writeup14b_10x_mouseembryo_dcca_laplacian_variablecalculations.R
