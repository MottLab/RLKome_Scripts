#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --account=def-gmott
#SBATCH --mem=124G
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=zachary.kileeg@mail.utoronto.ca
#SBATCH --mail-type=ALL


#load modules
module load StdEnv/2023
module load apptainer/1.2.4

#get variables from command line. Need the gff file containing features, the gff location for output, the genome
#being annotated, and the workdirectory for apptainer.
gff_in=${1} #name, not directory name
gff_out=${2}  #full name, not directory
genome_in=${3}
WORKDIR=${4}

#set variables for where to output the cds and peptides
cds_out="${gff_out}.cds.fasta"
prot_out="${gff_out}.pep.fasta"

#run apptainer. To get around weird file stuff, mount the workdirectory then add home to it. Apptainer will think
#the home directory is this newly mounted one and output stuff there.
apptainer run -C -W${SLURM_TMPDIR} -B${WORKDIR}:/home -B /scratch -B /localscratch /scratch/zkileeg/Ecotype_Data/AGAT/agat_1.2.0.sif \
	agat_sp_fix_overlaping_genes.pl --gff ${gff_in} -o "${gff_out}_merged.gff"  #merge overlapping loci
apptainer run -C -W${SLURM_TMPDIR} -B${WORKDIR}:/home -B /scratch -B /localscratch /scratch/zkileeg/Ecotype_Data/AGAT/agat_1.2.0.sif \
	agat_sp_keep_longest_isoform.pl -gff "${gff_out}_merged.gff" -o ${gff_out} #find representative sequence
apptainer run -C -W${SLURM_TMPDIR} -B${WORKDIR}:/home -B /scratch -B /localscratch /scratch/zkileeg/Ecotype_Data/AGAT/agat_1.2.0.sif \
	agat_sp_extract_sequences.pl -f $genome_in -g ${gff_out} -t cds -o ${cds_out} #output cds
apptainer run -C -W${SLURM_TMPDIR} -B${WORKDIR}:/home -B /scratch -B /localscratch /scratch/zkileeg/Ecotype_Data/AGAT/agat_1.2.0.sif \
	agat_sp_extract_sequences.pl -f $genome_in -g ${gff_out} -p -o ${prot_out}   #output peptides
