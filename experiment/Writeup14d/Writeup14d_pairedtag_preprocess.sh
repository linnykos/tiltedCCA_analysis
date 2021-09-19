#!/bin/bash
#$ -N pairedtag_preprocess
#$ -j y
#$ -o ../../../../out/Writeup14c/qsub/
#$ -l m_mem_free=300G

Rscript --no-save Writeup14d_pairedtag_preprocess.R
