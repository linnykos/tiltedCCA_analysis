#!/bin/bash
#$ -N citeseq_pbmc228
#$ -j y
#$ -o ../../../../out/Writeup14b/qsub/
#$ -l m_mem_free=200G

Rscript --no-save Writeup14b_citeseq_pbmc228_dcca.R
