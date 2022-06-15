#!/bin/bash
#$ -N abseq_97_distinct_de
#$ -j y
#$ -o ../../../../out/Writeup14n/qsub/
#$ -l m_mem_free=100G

Rscript --no-save abseq_bm_97antibodyRef_tcca_distinctDE.R
