#!/bin/bash
#SBATCH --time=0-48:00:00
#SBATCH --account=def-gmott
#SBATCH --mem=124G
#SBATCH --cpus-per-task=32
#SBATCH --mail-user=zachary.kileeg@mail.utoronto.ca
#SBATCH --mail-type=ALL


module load StdEnv/2020
module load samtools/1.16.1


sam_dir=${1}
out_dir=${2}


#Loop over the different arabidopsis ecotype files found within the accession file list.
for sam in $sam_dir*.sam; do

	if [ ! -d $out_dir ]; then
		mkdir $out_dir
	fi

       if [ ! -f "${out_dir}${sam##*/}.bam" ]; then

	#conert sam file to sorted bam file
		samtools view -@ 32 -bS -F 4 $sam | samtools sort -@ 32 -o "${out_dir}${sam##*/}.bam" -
	#index bam file
		samtools index -@ 32 "${out_dir}${sam##*/}.bam"
	fi

done
