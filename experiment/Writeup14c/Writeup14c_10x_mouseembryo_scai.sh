#!/bin/bash
#$ -N scai_mouseembryo
#$ -j y
#$ -o ../../../../out/Writeup14c/qsub/
#$ -l m_mem_free=200G

Rscript --no-save Writeup14c_10x_mouseembryo_scai.R
