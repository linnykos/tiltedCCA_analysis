#!/bin/bash
#$ -N citeseq_bm_gene
#$ -j y
#$ -o ../../../../out/Writeup14f/qsub/
#$ -l m_mem_free=50G

Rscript --no-save Writeup14f_citeseq_bm_dcca_postprocess_gene.R
