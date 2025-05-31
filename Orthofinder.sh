#!/bin/bash
#SBATCH --time=0-24:00:00
#SBATCH --account=def-gmott
#SBATCH --mem=124G
#SBATCH --cpus-per-task=32
#SBATCH --mail-user=zachary.kileeg@mail.utoronto.ca
#SBATCH --mail-type=ALL


module load StdEnv/2020 gcc/9.3.0
module load diamond/2.1.7
module load fastme/2.1.6.2
module load mcl/14.137


input_sequence_folder=${1}
output=${2}

/home/zkileeg/OrthoFinder/orthofinder -f ${input_sequence_folder} -o ${output} -y
