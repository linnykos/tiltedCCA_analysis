#!/bin/bash
#$ -N citeseq_pbmc224_dcca4
#$ -j y
#$ -o ../../../../out/Writeup14h/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup14h_citeseq_pbmc224_dcca4_zongming.R
