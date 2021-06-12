#!/bin/sh
#SBATCH -J 
#SBATCH -p
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --output
#SBATCH --error

export PATH=/data/mo/envs/Riboshape/bin:$PATH
echo ${PATH}
echo Begin quality control at `date`
sample=('sample_name_1' 'sample_name_2' 'sample_name_3')
path=('sample_path_1' 'sample_path_2' 'sample_path_3')

for i in $(seq 1 ${#sample[@]});do
echo Begin quality control of ${sample[$i-1]} at `date`
cd your_dir/
mkdir -p 0.fastqc/${sample[$i-1]}
fastqc \
your_dir/${sample[$i-1]}.clean.1.fastq.gz \
your_dir/${sample[$i-1]}.clean.2.fastq.gz \
--outdir your_dir/0.fastqc/${sample[$i-1]} --noextract >> your_dir/0.fastqc/${sample[$i-1]}/${sample[$i-1]}.log 2>&1
echo Done quality control of ${sample[$i-1]} at `date`
done

echo Done quanlity control at `date`
