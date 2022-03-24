#!/bin/bash
#$ -N greenleaf_tcca
#$ -j y
#$ -o ../../../../out/Writeup14m/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup14m_10x_greenleaf_activity_tcca.R
