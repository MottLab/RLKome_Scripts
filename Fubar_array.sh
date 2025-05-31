#!/bin/bash

#SBATCH --time=0-00:30:00
#SBATCH --account=def-gmott
#SBATCH --mem=24G
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=zachary.kileeg@mail.utoronto.ca
#SBATCH --mail-type=ALL
#SBATCH --array=1-591


module load StdEnv/2020  gcc/9.3.0  openmpi/4.0.3 hyphy/2.5.49

input_dir=${1}

mkdir FUBAR_output

trees=$(find $input_dir -maxdepth 1 -name "*.nwk" | sort > trees.txt )
alignments=$(find $input_dir -maxdepth 1 -name "*.fasta" | sort > alignments.txt)

tree_in=$( awk -v LineNum=$SLURM_ARRAY_TASK_ID 'NR==LineNum' trees.txt )
alignment=$( awk -v LineNum=$SLURM_ARRAY_TASK_ID 'NR==LineNum' alignments.txt )
#tree_in=$( awk -v LineNum=1 'NR==LineNum' trees.txt )
#alignment=$( awk -v LineNum=1 'NR==LineNum' alignments.txt )
echo $tree_in
echo $alignment

hyphy fubar --alignment ${alignment##*/}  --tree ${tree_in##*/} --output FUBAR_output/${alignment##*/}
