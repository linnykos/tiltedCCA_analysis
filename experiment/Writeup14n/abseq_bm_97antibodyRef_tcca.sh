#!/bin/bash
#$ -N abseq_97_tcca
#$ -j y
#$ -o ../../../../out/Writeup14n/qsub/
#$ -l m_mem_free=100G

Rscript --no-save abseq_bm_97antibodyRef_tcca.R
