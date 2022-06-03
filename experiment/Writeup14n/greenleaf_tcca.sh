#!/bin/bash
#$ -N greenleaf_tcca
#$ -j y
#$ -o ../../../../out/Writeup14n/qsub/
#$ -l m_mem_free=50G

Rscript --no-save greenleaf_tcca.R
