#!/bin/bash
#$ -N greenleaf_tcca
#$ -j y
#$ -o ../../../../out/Writeup14p/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup14p_greenleaf_tcca.R
