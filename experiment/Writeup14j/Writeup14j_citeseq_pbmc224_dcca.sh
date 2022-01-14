#!/bin/bash
#$ -N citeseq_pbmc224_dcca
#$ -j y
#$ -o ../../../../out/Writeup14j/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup14j_citeseq_pbmc224_dcca.R
