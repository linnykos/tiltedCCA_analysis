#!/bin/bash
#$ -N reik_tiltedcca3
#$ -j y
#$ -o ../../../../out/Writeup14p/qsub/
#$ -l m_mem_free=200G

Rscript --no-save Writeup14p_reik_tcca_v3.R
