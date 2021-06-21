#!/bin/sh
#SBATCH -J ribo_1_fastqc
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --output=ribo.1.fastq.out
#SBATCH --error=ribo.1.fastq.err

export PATH=your_envs_dir/:$PATH
echo start riboseq at `date`

sample=(sample1 sample2 sample3)
pth=()


echo start quality control at `date`
for i in $(seq 1 ${#sample[@]});do
mkdir -p your_output_dir/1.fastqc/${sample[$i-1]}
cd your_output_dir/1.fastqc/${sample[$i-1]}
fastqc -t 4 your_dir/${sample[$i-1]}.fq.gz --outdir ./ --noextract
echo finish ${sample[$i-1]} at `date`
done
echo finish quality control at `date`


echo start trim at `date`
for i in $(seq 1 ${#sample[@]});do
mkdir -p your_output_dir/2.trim/${sample[$i-1]}
cd your_output_dir/2.trim/${sample[$i-1]}
fastp --length_limit 50 \
--adapter_fasta your_dir/adapters.fa \
-i your_output_dir/1.fastqc/${sample[$i-1]}/${sample[$i-1]}.fq.gz \
-o your_output_dir/2.trim/${sample[$i-1]}/${sample[$i-1]}.clean.fq.gz \
--thread=4 -l 15 -j ${sample[$i-1]}.json -h ${sample[$i-1]}.html
echo finish trim ${sample[$i-1]} at `date`
done
echo finish trim at `date`


echo start remove rRNA at `date`
for i in $(seq 1 ${#sample[@]});do
echo start remove rRNA ${sample[$i-1]} at `date`
mkdir -p your_output_dir/3.remove_rRNA/fastq/${sample[$i-1]}
mkdir -p your_output_dir/3.remove_rRNA/rRNA/${sample[$i-1]}
bowtie -y -a --norc --best --strata -S -p 4 -l 15 \
--un your_output_dir/3.remove_rRNA/rRNA/${sample[$i-1]}/${sample[$i-1]}.rm_rRNA.fq \
your_dir/Arabidopsis_thaliana.TAIR10.34.rRNA \
-q your_output_dir/3.remove_rRNA/fastq/${sample[$i-1]}/${sample[$i-1]}.clean.fq.gz \
your_output_dir/3.remove_rRNA/fastq/${sample[$i-1]}/${sample[$i-1]}.alngned_rRNA.txt

cd your_output_dir/3.remove_rRNA/fastq/${sample[$i-1]}/
gzip ${sample[$i-1]}.rm_rRNA.fq
rm your_output_dir/3.remove_rRNA/fastq/${sample[$i-1]}/${sample[$i-1]}.alngned_rRNA.txt
echo finish remove rRNA ${sample[$i-1]} at `date`
done
echo finish remove rRNA at `date`
