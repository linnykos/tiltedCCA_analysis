#!/bin/bash
#$ -N greenleaf_timeseries
#$ -j y
#$ -o ../../../../out/Writeup14o/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup14o_10x_greenleaf_timeseries.R
