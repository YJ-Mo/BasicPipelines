#!/bin/sh
#SBATCH -J 1_fastqc_all
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --output=/home/user_25/scripts/log/1.fastqc.out
#SBATCH --error=/home/user_25/scripts/log/1.fastqc.err

export PATH=/data/zhaoyizi/software/anaconda3/envs/Riboshape/bin:$PATH
echo start qc `date`
sample=('CD1_1' 'CD1_2' 'CD1_3' 'CD0_1' 'CD0_2' 'CD0_3' 'CN1_1' 'CN1_2' 'CN1_3' 'CN0_1' 'CN0_2' 'CN0_3')
path=('WT/control/nouvb' 'WT/control/nouvb' 'WT/control/nouvb' 'WT/control/uvb' 'WT/control/uvb' 'WT/control/uvb' 'WT/modified/nouvb' 'WT/modified/nouvb' 'WT/modified/nouvb' 'WT/modified/uvb' 'WT/modified/uvb' 'WT/modified/uvb')
for i in $(seq 1 ${#sample[@]});do
cd /data/mo/PartI.RNA-seq
mkdir -p 1.fastqc/${sample[$i-1]}
fastqc \
/data/TA_QUIZ_RNA_regulation/data/riboshape_liulab_batch4/SHAPE-MaP/${path[$i-1]}/${sample[$i-1]}.clean.1.fastq.gz \
/data/TA_QUIZ_RNA_regulation/data/riboshape_liulab_batch4/SHAPE-MaP/${path[$i-1]}/${sample[$i-1]}.clean.2.fastq.gz \
--outdir /data/mo/PartI.RNA-seq/1.fastqc/${sample[$i-1]} --noextract >> /data/mo/PartI.RNA-seq/1.fastqc/${sample[$i-1]}/${sample[$i-1]}.log 2>&1
echo finish ${sample[$i-1]} `date`
done
echo finish qc `date`
