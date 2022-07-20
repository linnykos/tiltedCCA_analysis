#!/bin/bash
#$ -N reik_tiltedcca2
#$ -j y
#$ -o ../../../../out/Writeup14p/qsub/
#$ -l m_mem_free=200G

Rscript --no-save Writeup14p_reik_tcca_v2.R
