#!/bin/bash
#$ -N pbmc224_dcca_clisi
#$ -j y
#$ -o ../../../../out/Writeup14j/qsub/
#$ -l m_mem_free=200G

Rscript --no-save Writeup14j_citeseq_pbmc224_dcca_postprocess_clisi.R
