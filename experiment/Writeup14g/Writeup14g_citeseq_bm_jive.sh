#!/bin/bash
#$ -N citeseq_bm_jive
#$ -j y
#$ -o ../../../../out/Writeup14g/qsub/
#$ -l m_mem_free=50G

Rscript --no-save Writeup14g_citeseq_bm_jive.R
