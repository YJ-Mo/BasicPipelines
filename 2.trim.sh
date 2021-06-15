#!/bin/sh
#SBATCH -J 2_trim_CD0
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --output=/home/user_25/scripts/log/2.trimmed_CD0.out
#SBATCH --error=/home/user_25/scripts/log/2.trimmed_CD0.err

export PATH=/data/zhaoyizi/software/anaconda3/envs/Riboshape/bin:$PATH
echo start trimmed
sample=('CD0_1' 'CD0_2' 'CD0_3')
path=('WT/control/uvb' 'WT/control/uvb' 'WT/control/uvb')

for i in $(seq 1 ${#sample[@]});do
echo start ${sample[$i-1]} `date`
cd /data/mo/PartI.RNA-seq
mkdir -p 2.trim/${sample[$i-1]}
fastp -i /data/TA_QUIZ_RNA_regulation/data/riboshape_liulab_batch4/SHAPE-MaP/${path[$i-1]}/${sample[$i-1]}.clean.1.fastq.gz \
-I /data/TA_QUIZ_RNA_regulation/data/riboshape_liulab_batch4/SHAPE-MaP/${path[$i-1]}/${sample[$i-1]}.clean.2.fastq.gz \
-o /data/mo/PartI.RNA-seq/2.trim/${sample[$i-1]}/${sample[$i-1]}.clean.1.fastq.gz \
-O /data/mo/PartI.RNA-seq/2.trim/${sample[$i-1]}/${sample[$i-1]}.clean.2.fastq.gz \
--thread=4 -l 15
done
