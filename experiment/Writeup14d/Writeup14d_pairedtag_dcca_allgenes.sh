#!/bin/bash
#$ -N pairedtag_dcca_allgenes
#$ -j y
#$ -o ../../../../out/Writeup14d/qsub/
#$ -l m_mem_free=50G

Rscript --no-save Writeup14d_pairedtag_dcca_allgenes.R
