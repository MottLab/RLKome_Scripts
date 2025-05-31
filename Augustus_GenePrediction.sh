#!/bin/bash
#SBATCH --time=0-12:00:00
#SBATCH --account=def-gmott
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=zachary.kileeg@mail.utoronto.ca
#SBATCH --mail-type=ALL
#SBATCH --array=1-44


module load StdEnv/2020  gcc/9.3.0
module load augustus/3.4.0

genome_in_dir=${1}
output_dir=${2}

sleep $((SLURM_ARRAY_TASK_ID*10))


if [[ $SLURM_ARRAY_TASK_ID -eq 1 ]]; then

	find $genome_in_dir -type f -maxdepth 1 -name "*.fasta" | sort > genome_files.txt
fi

genome_in=$( awk -v LineNum=$SLURM_ARRAY_TASK_ID 'NR==LineNum' genome_files.txt )


output="${output_dir}${genome_in##*/}_augustus.gff"

echo "Query file ${genome_in}"
echo $output

#run augustus. settings are pretty much default, just set so the gff3 is the main output. It will be cleaned with the evidence modeler conversion script. 
augustus --gff3=on --UTR=on --outfile=$output --species=arabidopsis $genome_in

#cat $output | grep -Pzo "(?s)\[.*?\]" | \
#	sed -e 's/\#\s//g' | sed -e 's/\]\x00/\n/g' | awk '{gsub("\\[", ">Sequence\_"NR"\n", $0); print}' \
#	> $output".fasta"
