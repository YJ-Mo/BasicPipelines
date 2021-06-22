#!/bin/bash
#SBATCH -J 
#SBATCH -p 
#SBATCH --nodes=
#SBATCH --ntasks=
#SBATCH --output=your_dir/log.out
#SBATCH --error=your_dir/log.err

export PATH=your_envs_dir:$PATH

mkdir -p your_output_dir/
cd your_output_dir/

echo start calculating transcripts activity at `date`
your_software_dir/shapemap/shapemapper \
--target your_dir/ref.fa \.  # transcripts references fasta
--name sample_name \
--min-depth 100 \
--min-qual-to-count 20 \
--overwrite \
--modified --folder your_dir/  \ # path for folder containing data of modified condition
--untreated --folder your_dir/ \ # path for folder containing data of untreated condition
--star-aligner \
--nproc 8 \
--verbose

echo start calculating hit level at `date`
# hit level is calculated for each transcripts at different depth

path0=your_dir/ # dir containing hit_level.py
path1=your_dir/ # dir containing .hit files
path2=your_output_dir

## 注意，这里的输入文件是final.modified_unmodified而非final.modified_unmodified.hit
python $path0/hit_level.py \
--data_path $path1/final.modified_unmodified \
--savepath_hit $path2/final.modified_unmodified.hit
echo finish python script `date`

echo -e "cutoff\ttranscript_id\tmodified_depth_median\tunmodified_depth_median\tmodified_depth_sum\tunmodified_depth_sum\thit"       >       $path2/cutoff.hit.group;
awk -F '\t' '$3>0 && $2<=25 && $2>0{print "0\t"$0}'     $path2/final.modified_unmodified.hit  >>        $path2/cutoff.hit.group;
awk -F '\t' '$3>0 && $2<=50 && $2 > 25{print "25\t"$0}' $path2/final.modified_unmodified.hit  >>        $path2/cutoff.hit.group;
awk -F '\t' '$3>0 && $2<=100 && $2 > 50{print "50\t"$0}' $path2/final.modified_unmodified.hit >>      $path2/cutoff.hit.group;
awk -F '\t' '$3>0 && $2<=200 && $2 > 100{print "100\t"$0}' $path2/final.modified_unmodified.hit >>      $path2/cutoff.hit.group;
awk -F '\t' '$3>0 && $2<=300 && $2 > 200{print "200\t"$0}' $path2/final.modified_unmodified.hit >>      $path2/cutoff.hit.group;
awk -F '\t' '$3>0 && $2<=500 && $2 >300{print "300\t"$0}' $path2/final.modified_unmodified.hit >>      $path2/cutoff.hit.group;
awk -F '\t' '$3>0 && $2<=750 && $2 >500{print "500\t"$0}' $path2/final.modified_unmodified.hit >>      $path2/cutoff.hit.group;
awk -F '\t' '$3>0 && $2<=1000 && $2 >750{print "750\t"$0}' $path2/final.modified_unmodified.hit >>      $path2/cutoff.hit.group;
awk -F '\t' '$3>0 && $2<=2000 && $2 >1000{print "1000\t"$0}' $path2/final.modified_unmodified.hit >>      $path2/cutoff.hit.group;
awk -F '\t' '$3>0 && $2<=5000 && $2 >2000{print "2000\t"$0}' $path2/final.modified_unmodified.hit >>      $path2/cutoff.hit.group;
awk -F '\t' '$3>0 && $2 > 5000{print "5000\t"$0}' $path2/final.modified_unmodified.hit >>      $path2/cutoff.hit.group;

echo finish hit.group
