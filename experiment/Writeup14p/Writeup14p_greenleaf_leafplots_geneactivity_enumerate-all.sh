#!/bin/bash
#$ -N greenleaf_leafplot_enumerate
#$ -j y
#$ -o ../../../../out/Writeup14p/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup14p_greenleaf_leafplots_geneactivity_enumerate-all.R
