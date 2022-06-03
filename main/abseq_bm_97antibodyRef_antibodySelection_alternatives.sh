#!/bin/bash
#$ -N abseq_97ref_varSelect_alternative
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=100G

Rscript --no-save abseq_bm_97antibodyRef_antibodySelection_alternatives.R
