#!/bin/bash
#$ -N dcca_H3K4me1
#$ -j y
#$ -o ../../../../out/Writeup14c/qsub/
#$ -l m_mem_free=300G

Rscript --no-save Writeup14c_pairedtag_H3K4me1_dcca.R
