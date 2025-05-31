#!/bin/bash

#SBATCH --time=2-00:00:00
#SBATCH --account=def-gmott
#SBATCH --mem=124G
#SBATCH --cpus-per-task=32
#SBATCH --mail-user=zachary.kileeg@mail.utoronto.ca
#SBATCH --mail-type=ALL
#SBATCH --array=2-2

module load StdEnv/2023
module load hisat2/2.2.1


paired_end_fastq_dir=${1}
single_end_fastq_dir=${2}
genome_folder_in=${3}
mapped_output_dir=${4}



#make temp directory to hold sam files
mkdir "${mapped_output_dir}temp_${SLURM_ARRAY_TASK_ID}"
cd "${mapped_output_dir}temp_${SLURM_ARRAY_TASK_ID}"

#get comma separated data for each fastq file in the paired end and single end  RNA seq folders
fastq_PE_1=$( find $paired_end_fastq_dir -type f -maxdepth 1 -name "*_1.fastq" | paste -s -d, - )
fastq_PE_2=$( find $paired_end_fastq_dir -type f -maxdepth 1 -name "*_2.fastq" | paste -s -d, - )
fastq_SE=$( find $single_end_fastq_dir -type f -maxdepth 1 -name "*.fastq" | paste -s -d, - )

if [ ! -f genome_files.txt ]; then
	find $genome_folder_in -type f -name "*.fasta" -maxdepth 1 | sort > genome_files.txt
fi
genome_file_in=$( awk -v LineNum=$SLURM_ARRAY_TASK_ID 'NR==LineNum' genome_files.txt )

#echo $paired_end_fastq_dir
#echo $single_end_fastq_dir
#echo $genome_file_in
#echo $mapped_output_dir

#echo "Test as string variable = ${fastq_array_find1}"
#echo "Test as array = ${fastq_array_find1[@]}"

hisat2-build -p 32 $genome_file_in "${genome_file_in##*/}_index"              #build index with the genome file, outputting

hisat2 -p 32 -x "${genome_file_in##*/}_index" -U $fastq_SE -S "${genome_file_in##*/}_SE_Hisat2.sam"

samtools view -@ 32 -bS -F 4 "${genome_file_in##*/}_SE_Hisat2.sam" | \
samtools sort -@ 32 -o "${mapped_output_dir}${genome_file_in##*/}_SE_Hisat2.bam" -
#index bam file
samtools index -@ 32 "${mapped_output_dir}${genome_file_in##*/}_SE_Hisat2.bam"

cd ${mapped_output_dir}
rm -r temp_${SLURM_ARRAY_TASK_ID}


