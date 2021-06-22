#!/bin/bash
#SBATCH -J 
#SBATCH -p 
#SBATCH --nodes=
#SBATCH --ntasks=
#SBATCH --output=your_dir/log.out
#SBATCH --error=your_dir/log.err

export PATH=your_envs_dir:$PATH

echo start calculating transcripts activity
your_software_dir/shapemap/shapemapper \
--target your_dir/ref.fa \.  # transcripts references fasta
--name sample_name \
--min-depth 100 \
--min-qual-to-count 20 \
--overwrite \
--modified --folder your_dir/  \
--untreated --folder your_dir/ \
--star-aligner \
--nproc 8 \
--verbose

echo
