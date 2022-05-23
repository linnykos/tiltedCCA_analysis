#!/bin/bash
#$ -N greenleaf_enumerate
#$ -j y
#$ -o ../../../../out/Writeup14n/qsub/
#$ -l m_mem_free=50G

Rscript --no-save Writeup14n_greenleaf_enumerate.R
