#!/bin/bash
#$ -N reik_enumerate
#$ -j y
#$ -o ../../../../out/Writeup14p/qsub/
#$ -l m_mem_free=200G

Rscript --no-save Writeup14p_reik_enumerate.R
