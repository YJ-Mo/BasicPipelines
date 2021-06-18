#!/bin/sh
#SBATCH -J #task name
#SBATCH -p #node name
#SBATCH --nodes= #num of node
#SBATCH --ntasks= #num of task
#SBATCH --output=your_dir/your_task.out
#SBATCH --error=your_dir/your_task.err

export PATH=your_dir/:$PATH #path for envs
samples=(your_sample)
paths=(your_sample_path)

echo start quality control at `date`
for i in $(seq 1 ${#samples[@]});do
mkdir -p your_output_dir/0.fastqc/${sample[$i-1]}
cd your_output_dir/1.fastqc/${sample[$i-1]}
fastqc \
your_dir/${path[$i-1]}/${sample[$i-1]}.clean.1.fastq.gz \
your_dir/${path[$i-1]}/${sample[$i-1]}.clean.2.fastq.gz \
--outdir ./ --noextract >> your_dir/0.fastqc/${sample[$i-1]}/${sample[$i-1]}.log 2>&1
echo finish QC of ${sample[$i-1]} at `date`
done
echo finish all qc at `date`

echo start trimming at `date`
for i in $(seq 1 ${#samples[@]});do
mkdir -p your_output_dir/2.trim/${sample[$i-1]}
cd your_output_dir/2.trim/${sample[$i-1]}
fastp -i your_dir/${path[$i-1]}/${sample[$i-1]}.clean.1.fastq.gz \
-I your_dir/${path[$i-1]}/${sample[$i-1]}.clean.2.fastq.gz \
-o ./${sample[$i-1]}.clean.1.fastq.gz \
-O ./${sample[$i-1]}.clean.2.fastq.gz \
--thread=4 -l 15
echo finish trimming ${sample[$i-1]} at `date`
done
echo finish all trimming at `date`

echo start remove rRNA at `date`
for i in $(seq 1 ${#samples[@]});do
mkdir -p your_output_dir/3.remove_rRNA/fastq/${sample[$i-1]}
mkdir -p your_output_dir/3.remove_rRNA/rRNA/${sample[$i-1]}
bowtie -n 0 -y -a --norc --best --strata -S -p 4 -l 15 \
--un your_output_dir/3.remove_rRNA/fastq/${sample[$i-1]}/${sample[$i-1]}.rm_rRNA.fq \
your_dir/Arabidopsis_thaliana.TAIR10.34.rRNA \
-1 your_output_dir/2.trim/${sample[$i-1]}/${sample[$i-1]}.clean.1.fastq.gz \
-2 your_output_dir/2.trim/${sample[$i-1]}/${sample[$i-1]}.clean.2.fastq.gz \
your_output_dir/3.remove_rRNA/rRNA/${sample[$i-1]}/${sample[$i-1]}.aligned_rRNA.txt
echo finish remove rRNA ${sample[$i-1]} at `date`

rm your_output_dir/2.trim/${sample[$i-1]}/${sample[$i-1]}.clean.1.fastq.gz
rm your_output_dir/2.trim/${sample[$i-1]}/${sample[$i-1]}.clean.2.fastq.gz
echo rm trim files at `date`

cd your_output_dir/3.remove_rRNA/fastq/${sample[$i-1]}/
gzip ${sample[$i-1]}.rm_rRNA_1.fq
gzip ${sample[$i-1]}.rm_rRNA_2.fq
rm ${sample[$i-1]}.rm_rRNA_1.fq
rm ${sample[$i-1]}.rm_rRNA_2.fq
echo finish gzip ${sample[$i-1]} at `date`
echo finish remove rRNA ${sample[$i-1]} `date`
done
echo finish all remove rRNA at `date`

echo start mapping at `date`
for i in $(seq 1 ${#samples[@]});do
mkdir -p your_output_dir/4.mapping/${sample[$i-1]}
cd your_output_dir/4.mapping/${sample[$i-1]}
STAR \
--runThreadN 8 \
--limitBAMsortRAM 60000000000 \
--outFilterType BySJout \
--outFilterMismatchNmax 10  \
--genomeDir your_dir/genome \
--readFilesIn your_output_dir//3.remove_rRNA/fastq/${sample[$i-1]}/${sample[$i-1]}.rm_rRNA_1.fq.gz \
your_output_dir//3.remove_rRNA/fastq/${sample[$i-1]}/${sample[$i-1]}.rm_rRNA_2.fq.gz \
--readFilesCommand 'zcat' \
--outFileNamePrefix  ${sample[$i-1]} \
--outSAMtype BAM Unsorted \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMattributes All  --outSAMstrandField intronMotif --outBAMcompression 6 --outReadsUnmapped Fastx
echo finish mapping ${sample[$i-1]} at `date`

samtools sort -T \
your_output_dir/4.mapping/${sample[$i-1]}/${sample[$i-1]}_Aligned.out.sorted \
-o your_output_dir/4.mapping/${sample[$i-1]}/${sample[$i-1]}_Aligned.sortedByCoord.out.bam \
your_output_dir/4.mapping/${sample[$i-1]}/${sample[$i-1]}_Aligned.out.bam

samtools sort -T \
your_output_dir/4.mapping/${sample[$i-1]}/${sample[$i-1]}_Aligned.toTranscriptome.out.sorted \
-o your_output_dir/4.mapping/${sample[$i-1]}/${sample[$i-1]}_Aligned.toTranscriptome.out.sorted.bam \
your_output_dir/4.mapping/${sample[$i-1]}/${sample[$i-1]}_Aligned.toTranscriptome.out.bam

samtools index \
your_output_dir/4.mapping/${sample[$i-1]}/${sample[$i-1]}_Aligned.sortedByCoord.out.bam

samtools index \
your_output_dir/4.mapping/${sample[$i-1]}/${sample[$i-1]}_Aligned.toTranscriptome.out.sorted.bam

echo finish samtool ${sample[$i-1]} at `date`
done
echo finish all mapping at `date`

echo start counting at `date`
for i in $(seq 1 ${#samples[@]});do
mkdir -p your_output_dir/5.read_counts/${sample[$i-1]}
cd your_output_dir/5.read_counts/${sample[$i-1]}

featureCounts \
-T 8 \
-s 0 \
-p -t CDS \
-g gene_id \
-a your_dir/Arabidopsis_thaliana.TAIR10.34.gtf \
-o your_output_dir/5.read_counts/${sample[$i-1]}/${sample[$i-1]}_featurecounts.txt \
your_output_dir/4.mapping/${sample[$i-1]}/${sample[$i-1]}_Aligned.sortedByCoord.out.bam

featureCounts \
-T 8 \
-s 0 \
-p -t exon \
-g gene_id \
-a your_dir/Arabidopsis_thaliana.TAIR10.34.gtf \
-o your_output_dir/5.read_counts/${sample[$i-1]}/${sample[$i-1]}_featurecounts.all.txt \
your_output_dir/4.mapping/${sample[$i-1]}/${sample[$i-1]}_Aligned.sortedByCoord.out.bam

mkdir -p your_output_dir/5.read_counts/result/${sample[$i-1]}
cd your_output_dir/5.read_counts/result/${sample[$i-1]}
echo -e "gene_id	${sample[$i-1]}" >your_output_dir/5.read_counts/result/${sample[$i-1]}.txt
cat your_output_dir/5.read_counts/${sample[$i-1]}/${sample[$i-1]}_featurecounts.txt| grep -v '#' | grep -v 'Geneid' | cut -f 1,7 >> your_output_dir/5.read_counts/result/${sample[$i-1]}.txt

echo -e "gene_id	${sample[$i-1]}" >your_output_dir/5.read_counts/result/${sample[$i-1]}.all.txt
cat your_output_dir/5.read_counts/${sample[$i-1]}/${sample[$i-1]}_featurecounts.all.txt| grep -v '#' | grep -v 'Geneid' | cut -f 1,7 >> your_output_dir/5.read_counts/result/${sample[$i-1]}.all.txt

echo finish counting ${sample[$i-1]} at `date`
done 
echo finish all counting at `date`

echo start merge counting at `date`
for i in $(seq 1 ${#samples[@]});do
mkdir -p your_output_dir/6.mergeCount
python your_dir/mergeCount.py
echo finish merge counting ${sample[$i-1]} at `date`
done
echo finish merge counting at `date`


