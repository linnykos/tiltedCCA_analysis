#!/bin/bash
#$ -N citeseq_pbmc25
#$ -j y
#$ -o ../../../../out/Writeup14b/qsub/
#$ -l m_mem_free=50G

Rscript --no-save Writeup14b_citeseq_pbmc25_dcca.R
