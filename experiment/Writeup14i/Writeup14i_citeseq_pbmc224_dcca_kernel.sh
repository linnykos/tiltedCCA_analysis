#!/bin/bash
#$ -N citeseq_pbmc224_dcca
#$ -j y
#$ -o ../../../../out/Writeup14i/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup14i_citeseq_pbmc224_dcca_kernel.R