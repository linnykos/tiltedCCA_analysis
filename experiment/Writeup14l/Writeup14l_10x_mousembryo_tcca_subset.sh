#!/bin/bash
#$ -N mouse_10x_tcca_subset
#$ -j y
#$ -o ../../../../out/Writeup14l/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup14l_10x_mousembryo_tcca_subset.R
