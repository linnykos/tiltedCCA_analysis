#!/bin/bash
#$ -N pairedtag_dcca
#$ -j y
#$ -o ../../../../out/Writeup14d/qsub/
#$ -l m_mem_free=50G

Rscript --no-save Writeup14d_pairedtag_dcca.R
