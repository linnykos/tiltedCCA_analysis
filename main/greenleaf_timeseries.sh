#!/bin/bash
#$ -N greenleaf_timeseries
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=100G

Rscript --no-save greenleaf_timeseries.R
