#!/bin/bash
#$ -N greenleaf_10x_tcca
#$ -j y
#$ -o ../../../../out/Writeup14l/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup14l_10x_greenleaf_tcca.R
