#!/bin/bash
#$ -N SNU601_atac
#$ -j y
#$ -o ../../../../out/Writeup14e/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup14e_SNU601_atac_exploratory.R
