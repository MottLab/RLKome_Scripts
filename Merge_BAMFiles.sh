#!/bin/bash
#SBATCH --time=0-24:00:00
#SBATCH --account=def-gmott
#SBATCH --mem=124G
#SBATCH --cpus-per-task=32
#SBATCH --mail-user=zachary.kileeg@mail.utoronto.ca
#SBATCH --mail-type=ALL


module load StdEnv/2020
module load samtools/1.16.1


bam_dir=${1}
output_name=${2}

#Loop over the different arabidopsis ecotype files found within the accession file list.


bam_files=$( find ${bam_dir} -type f  -name "*.bam" )
#echo ${bam_files[@]}
#conert sam file to sorted bam file
samtools merge -@ 32 -f -o ${output_name} $bam_files
#index bam file
samtools index -@ 32 ${output_name}



