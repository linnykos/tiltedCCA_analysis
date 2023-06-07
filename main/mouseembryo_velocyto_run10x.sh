#!/bin/bash/
#$ -N velo_run10x_mouseembryo
#$ -j y
#$ -o /home/stat/kevinl1/project/tiltedCCA/out/main/qsub
#$ -l m_mem_free=100G

source $HOME/project/Multiome_fate/venv396/bin/activate
module load samtools

velocyto run10x -m $HOME/nzhanglab/data/reference/repeat_annotation/mm10_rmsk.gtf $HOME/nzhanglab/data/10x_mouse_embryo $HOME/nzhanglab/data/reference/refdata-gex-mm10-2020-A/genes/genes.gtf
