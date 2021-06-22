#!/bin/bash
#SBATCH -J 
#SBATCH -p 
#SBATCH --nodes=
#SBATCH --ntasks=
#SBATCH --output=/shape.0.map.out
#SBATCH --error=/shape.0.map.err

export PATH=/data/zhaoyizi/software/anaconda3/envs/Riboshape/bin:/data/liuxiaofan/software/shapemap/:$PATH
cd /data/mo/SHAPE

/data/liuxiaofan/software/shapemap/shapemapper \
--target /data/TA_QUIZ_RNA_regulation/data/test/Arabidopsis_thaliana.TAIR10.34.transcripts_new_2.fa \
--name "C1_1" \
--min-depth 100 \
--min-qual-to-count 20 \
--overwrite \
--modified --folder /data/TA_QUIZ_RNA_regulation/data/test/modified  \
--untreated --folder /data/TA_QUIZ_RNA_regulation/data/test/control \
--star-aligner \
--nproc 8 \
--verbose
