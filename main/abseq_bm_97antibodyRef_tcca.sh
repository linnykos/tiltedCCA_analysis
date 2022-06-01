#!/bin/bash
#$ -N abseq_97ref_tcca
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=50G

Rscript --no-save abseq_bm_97antibodyRef_tcca.R
