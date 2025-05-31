#!/bin/bash

#SBATCH --time=0-24:00:00
#SBATCH --account=def-gmott
#SBATCH --mem=124G
#SBATCH --cpus-per-task=32
#SBATCH --mail-user=zachary.kileeg@mail.utoronto.ca
#SBATCH --mail-type=ALL
#SBATCH --array=1-44

module load StdEnv/2023
module load hisat2/2.2.1


fasta_reference_folder=${1}
tree_folder=${2}
output_dir=${3}

find $fasta_reference_folder -type f -name "*.fasta" | sort > genome_files.txt

fasta_reference=$( awk -v LineNum=$SLURM_ARRAY_TASK_ID 'NR==LineNum' genome_files.txt )
fasta_index_name="${fasta_reference##*/}"

find $tree_folder -type f -name "*.nwk" | sort > genome_files.txt

tree_in=$( awk -v LineNum=$SLURM_ARRAY_TASK_ID 'NR==LineNum' genome_files.txt )

hyphy FUBAR --alignment $fasta_reference --tree $tree_in --output ${output_dir}
