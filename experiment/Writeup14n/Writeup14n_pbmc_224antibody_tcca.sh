#!/bin/bash
#$ -N pbmc_224
#$ -j y
#$ -o ../../../../out/Writeup14n/qsub/
#$ -l m_mem_free=200G

Rscript --no-save Writeup14n_pbmc_224antibody_tcca.R
