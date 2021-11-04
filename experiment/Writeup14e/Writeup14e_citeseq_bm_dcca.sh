#!/bin/bash
#$ -N citeseq
#$ -j y
#$ -o ../../../../out/Writeup14e/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup14e_citeseq_bm_dcca.R
