#!/bin/bash
#$ -N dcca_H3K9me3
#$ -j y
#$ -o ../../../../out/Writeup14c/qsub/
#$ -l m_mem_free=200G

Rscript --no-save Writeup14c_pairedtag_H3K9me3_dcca.R
