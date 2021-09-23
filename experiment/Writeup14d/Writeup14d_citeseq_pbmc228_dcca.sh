#!/bin/bash
#$ -N citeseq_pbmc_dcca
#$ -j y
#$ -o ../../../../out/Writeup14d/qsub/
#$ -l m_mem_free=200G

Rscript --no-save Writeup14d_citeseq_pbmc228_dcca.R
