#!/bin/bash
#$ -N pbmc224_dcca2
#$ -j y
#$ -o ../../../../out/Writeup14j/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup14j_citeseq_pbmc224_dcca2.R
