#!/bin/sh
#SBATCH -J mapping
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --output=/home/user_25/scripts/log/5.mapping_CD0.out
#SBATCH --error=/home/user_25/scripts/log/5.mapping_CD0.err

export PATH=/data/mo/envs:$PATH
sample=('CD0_1' 'CD0_2' 'CD0_3')

for i in $(seq 1 ${#sample[@]});do
echo mapping ${sample[$i-1]} `date`
cd /data/mo/PartI.RNA-seq
mkdir -p 4.mapping/${sample[$i-1]}
cd 4.mapping/${sample[$i-1]}

STAR \
--runThreadN 8 \
--limitBAMsortRAM 60000000000 \
--outFilterType BySJout \
--outFilterMismatchNmax 10  \
--genomeDir /data/TA_QUIZ_RNA_regulation/data/ATH/index/STAR/genome \
--readFilesIn /data/mo/PartI.RNA-seq/3.remove_rRNA/fastq/${sample[$i-1]}/${sample[$i-1]}.rm_rRNA_1.fq.gz \
/data/mo/PartI.RNA-seq/3.remove_rRNA/fastq/${sample[$i-1]}/${sample[$i-1]}.rm_rRNA_2.fq.gz \
--readFilesCommand 'zcat' \
--outFileNamePrefix  ${sample[$i-1]} \
--outSAMtype BAM Unsorted \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMattributes All  --outSAMstrandField intronMotif --outBAMcompression 6 --outReadsUnmapped Fastx

echo finish mapping ${sample[$i-1]} `date`

samtools sort -T \
/data/mo/PartI.RNA-seq/4.mapping/${sample[$i-1]}/${sample[$i-1]}Aligned.out.sorted \
-o /data/mo/PartI.RNA-seq/4.mapping/${sample[$i-1]}/${sample[$i-1]}Aligned.sortedByCoord.out.bam \
/data/mo/PartI.RNA-seq/4.mapping/${sample[$i-1]}/${sample[$i-1]}Aligned.out.bam

samtools sort -T \
/data/mo/PartI.RNA-seq/4.mapping/${sample[$i-1]}/${sample[$i-1]}Aligned.toTranscriptome.out.sorted \
-o /data/mo/PartI.RNA-seq/4.mapping/${sample[$i-1]}/${sample[$i-1]}Aligned.toTranscriptome.out.sorted.bam \
/data/mo/PartI.RNA-seq/4.mapping/${sample[$i-1]}/${sample[$i-1]}Aligned.toTranscriptome.out.bam

samtools index \
/data/mo/PartI.RNA-seq/4.mapping/${sample[$i-1]}/${sample[$i-1]}Aligned.sortedByCoord.out.bam

samtools index \
/data/mo/PartI.RNA-seq/4.mapping/${sample[$i-1]}/${sample[$i-1]}Aligned.toTranscriptome.out.sorted.bam

echo finish samtool ${sample[$i-1]} `date`
done

echo mapping success `date`
