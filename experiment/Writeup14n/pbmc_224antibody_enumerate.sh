#!/bin/bash
#$ -N pbmc224_enumerate
#$ -j y
#$ -o ../../../../out/Writeup14n/qsub/
#$ -l m_mem_free=100G

Rscript --no-save pbmc_224antibody_enumerate.R
