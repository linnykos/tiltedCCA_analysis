#!/bin/bash
#$ -N H3K9me3_dcca
#$ -j y
#$ -o ../../../../out/Writeup14b/qsub/
#$ -l m_mem_free=50G

Rscript --no-save Writeup14b_pairedtag_H3K9me3_dcca.R
