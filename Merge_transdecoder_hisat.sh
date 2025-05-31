#!/bin/bash
#SBATCH --time=0-00:10:00
#SBATCH --account=def-gmott
#SBATCH --mem=24G
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=zachary.kileeg@mail.utoronto.ca
#SBATCH --mail-type=ALL




transdecoder_suffix=${1}
hisat_suffix=${2}
folder=${3}
output_suffix=${4}

cat ${folder}${folder%*/}${transdecoder_suffix} ${folder}${folder%*/}${hisat_suffix} > ${folder}${folder%*/}${output_suffix}



