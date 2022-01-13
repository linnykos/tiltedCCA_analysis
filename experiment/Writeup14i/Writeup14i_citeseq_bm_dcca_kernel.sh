#!/bin/bash
#$ -N citeseq_bm_dcca
#$ -j y
#$ -o ../../../../out/Writeup14i/qsub/
#$ -l m_mem_free=50G

Rscript --no-save Writeup14i_citeseq_bm_dcca_kernel.R
