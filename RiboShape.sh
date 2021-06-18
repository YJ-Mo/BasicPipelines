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
