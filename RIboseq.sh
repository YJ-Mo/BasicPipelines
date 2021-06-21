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
