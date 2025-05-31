#!/bin/bash
#SBATCH --time=0-12:00:00
#SBATCH --account=def-gmott
#SBATCH --mem=164G
#SBATCH --cpus-per-task=32
#SBATCH --mail-user=zachary.kileeg@mail.utoronto.ca
#SBATCH --mail-type=ALL
#SBATCH --array=1-5


#Exonerate is the bane of my existance. To whoever is reading this...don't use exonerate 2.4.0. The multi-core option is not stable and it will crash 
#a lot. I think it's because of the cpu architecture but not sure. 
module load StdEnv/2020  gcc/9.3.0
module load exonerate/2.4.0

#get input from cmd line
protein_reference=${1}
target_assembly_folder=${2}
output_folder=${3}

if [ "${output_folder: -1}" != "/" ]; then
	output_folder="${output_folder}/"
fi

find ${target_assembly_folder} -type f -name *.fasta -maxdepth 1 | sort > genome_files.txt

target_assembly=$(awk -v LineNum=$SLURM_ARRAY_TASK_ID 'NR==LineNum' genome_files.txt )

echo $protein_reference
echo $target_assembly
echo $output_folder

output_name="${target_assembly##*/}_exonerate.gff"
#run exonerate. homologous match must have 70 percent identity, introns must be between 10-60000 bp. Output as gff. run in mode protein2genome. 
exonerate -m protein2genome -c 32 --querytype protein --targettype dna --percent 70 --minintron 10 --maxintron 60000 --showtargetgff TRUE --showvulgar no --showalignment no \
 $protein_reference $target_assembly > "${output_folder}${output_name}"
