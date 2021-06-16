#!/bin/sh
#SBATCH -J 3_removerRNA_CD0
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --output=/home/user_25/scripts/log/3.rmrRNA_CD0.out
#SBATCH --error=/home/user_25/scripts/log/3.rmrRNA_CD0.err

export PATH=/data/mo/env/bin:$PATH

sample=('CD0_1' 'CD0_2' 'CD0_3')
path=('WT/control/uvb' 'WT/control/uvb' 'WT/control/uvb')

echo start CD0 `date`

for i in $(seq 1 ${#sample[@]});do
echo start ${sample[$i-1]} `date`
cd /data/mo/PartI.RNA-seq
mkdir -p 3.remove_rRNA/fastq/${sample[$i-1]}
mkdir -p 3.remove_rRNA/rRNA/${sample[$i-1]}
bowtie -n 0 -y -a --norc --best --strata -S -p 4 -l 15 \
--un /data/mo/PartI.RNA-seq/3.remove_rRNA/fastq/${sample[$i-1]}/${sample[$i-1]}.rm_rRNA.fq \
/data/TA_QUIZ_RNA_regulation/data/ATH/index/bowtie1/rRNA/Arabidopsis_thaliana.TAIR10.34.rRNA \
-1 /data/mo/PartI.RNA-seq/2.trim/${sample[$i-1]}/${sample[$i-1]}.clean.1.fastq.gz \
-2 /data/mo/PartI.RNA-seq/2.trim/${sample[$i-1]}/${sample[$i-1]}.clean.2.fastq.gz \
/data/mo/PartI.RNA-seq/3.remove_rRNA/rRNA/${sample[$i-1]}/${sample[$i-1]}.aligned_rRNA.txt
echo finish bowtie
rm /data/mo/PartI.RNA-seq/2.trim/${sample[$i-1]}/${sample[$i-1]}.clean.1.fastq.gz
rm /data/mo/PartI.RNA-seq/2.trim/${sample[$i-1]}/${sample[$i-1]}.clean.2.fastq.gz
echo finish rm

cd /data/mo/PartI.RNA-seq/3.remove_rRNA/fastq/${sample[$i-1]}/
gzip ${sample[$i-1]}.rm_rRNA_1.fq
gzip ${sample[$i-1]}.rm_rRNA_2.fq
echo finish gzip
rm ${sample[$i-1]}.rm_rRNA_1.fq
rm ${sample[$i-1]}.rm_rRNA_2.fq
echo finish rm
echo finish ${sample[$i-1]} `date`
done
echo remove_rRNA success `date`
