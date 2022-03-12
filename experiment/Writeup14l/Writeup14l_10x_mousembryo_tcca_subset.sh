#!/bin/bash
#$ -N 10x_mouse_tcca_subset
#$ -j y
#$ -o ../../../../out/Writeup14l/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup14l_10x_mousembryo_tcca_subset.R
