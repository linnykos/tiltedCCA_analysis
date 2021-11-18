#!/bin/bash
#$ -N SNU601_dcca_clisi
#$ -j y
#$ -o ../../../../out/Writeup14f/qsub/
#$ -l m_mem_free=50G

Rscript --no-save Writeup14f_SNU601_dcca_postprocess_clisi.R
