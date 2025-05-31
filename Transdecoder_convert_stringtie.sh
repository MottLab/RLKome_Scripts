#!/bin/bash
#SBATCH --time=0-1:00:00
#SBATCH --account=def-gmott
#SBATCH --mem=124G
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=zachary.kileeg@mail.utoronto.ca
#SBATCH --mail-type=ALL

module load StdEnv/2020  gcc/9.3.0
module load transdecoder/5.7.1

input_gtf=${1}
input_genome=${2}
output_transcripts=${3}
library_name=${4}

output_transcripts_transdecoder="${output_transcripts}.transdecoder.gff3"
stringtie_gff="${input_gtf}.gff"

#take input as gtf file from stringtie, input fasta genome, and output trancript prediction as fasta
gtf_genome_to_cdna_fasta.pl $input_gtf $input_genome > $output_transcripts

#take as input gtf file from stringtie and output as gff
gtf_to_alignment_gff3.pl $input_gtf > $stringtie_gff

#find ORFs using transcripts in fasta format and output as a gff3 file
TransDecoder.LongOrfs -t $output_transcripts
TransDecoder.Predict -t $output_transcripts

#get output predicted as gff3 file for the genome
cdna_alignment_orf_to_genome_orf.pl "${library_name}_transcripts.fasta.transdecoder_dir/longest_orfs.cds.best_candidates.gff3" \
	$stringtie_gff $output_transcripts > "${library_name}.genome.gff3"
