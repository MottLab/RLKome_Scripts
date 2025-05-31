#!/bin/bash
#SBATCH --time=0-02:00:00
#SBATCH --account=def-gmott
#SBATCH --mem=124G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-44
#SBATCH --mail-user=zachary.kileeg@mail.utoronto.ca
#SBATCH --mail-type=ALL


#load modules
module load StdEnv/2023
module load apptainer/1.2.4

#get variables from command line. Need the gff file containing features, the gff location for output, the genome
#being annotated, and the workdirectory for apptainer.
gff_folder_in=${1}
genome_folder=${2}
WORKDIR=${3}

#get a list of genome files to use and output it in the current directory 
find $genome_folder -type f -maxdepth 1 -name "*.fasta" | sort > genome_files.txt
find $gff_folder_in -type f -maxdepth 1 \( -iname \*.gff -o -iname \*_polished \)| sort > gff_files.txt
#this is a fun line that will get the value in the line using the internal slurm array task number as the line number
genome_in=$(awk -v LineNum=$SLURM_ARRAY_TASK_ID 'NR==LineNum' genome_files.txt)
#line=$(awk -v LineNum=1 'NR==LineNum' genome_files.txt)
gff_in=$(awk -v LineNum=$SLURM_ARRAY_TASK_ID 'NR==LineNum' gff_files.txt)

#set variables for where to output the cds and peptides

#substitution in variable to make life easier for pulling out only file name
#filename=$(basename "${line}" | sed 's/fasta/gff/g' )


#gff_in=${filename}
gff_merged="${WORKDIR}${gff_in##*/}_merged.gff"
gff_rep="${WORKDIR}${gff_in##*/}_representative.gff"
cds_out="${WORKDIR}${gff_in##*/}_representative.cds.fasta"
prot_out="${WORKDIR}${gff_in##*/}_representative.pep.fasta"

echo $filename
echo $genome_folder
echo $genome_in
echo $gff_in
echo $gff_merged
echo $cds_out
echo $prot_out
#run apptainer. To get around weird file stuff, mount the workdirectory then add home to it. Apptainer will think
#the home directory is this newly mounted one and output stuff there.
apptainer run -C -W${SLURM_TMPDIR} -B /home -B /scratch -B /localscratch /scratch/zkileeg/Ecotype_Data/AGAT/agat_1.2.0.sif \
	agat_sp_fix_overlaping_genes.pl --gff ${gff_in} -o ${gff_merged}  #merge overlapping loci
apptainer run -C -W${SLURM_TMPDIR} -B /home -B /scratch -B /localscratch /scratch/zkileeg/Ecotype_Data/AGAT/agat_1.2.0.sif \
	agat_sp_keep_longest_isoform.pl -gff ${gff_merged} -o ${gff_rep} #find representative sequence
apptainer run -C -W${SLURM_TMPDIR} -B /home -B /scratch -B /localscratch /scratch/zkileeg/Ecotype_Data/AGAT/agat_1.2.0.sif \
	agat_sp_extract_sequences.pl -f ${genome_in} -g ${gff_rep} -t cds -o ${cds_out} #output cds
apptainer run -C -W${SLURM_TMPDIR} -B /home -B /scratch -B /localscratch /scratch/zkileeg/Ecotype_Data/AGAT/agat_1.2.0.sif \
	agat_sp_extract_sequences.pl -f ${genome_in} -g ${gff_rep} -p -o ${prot_out}   #output peptides



#gff_in="${gff_folder_in}${line}.corr.Athaliana_447_Araport11.liftoff.gff3_polished"
#gff_merged="${gff_in}_merged.gff"
#cds_out="${gff_in}.cds.fasta"
#prot_out="${gff_in}.pep.fasta"


#gff_in="${gff_folder_in}${line}.corr_Athaliana_447_Araport11.liftoff.gff3_polished"
#gff_merged="${gff_in}_merged.gff"
#gff_rep="${gff_in}_representative.gff"
#cds_out="${gff_in}_representative.cds.fasta"
#prot_out="${gff_in}_representative.pep.fasta"

#echo $genome_in
#echo $gff_in
#echo $gff_merged
#echo $cds_out
#echo $prot_out
#run apptainer. To get around weird file stuff, mount the workdirectory then add home to it. Apptainer will think
#the home directory is this newly mounted one and output stuff there.
#apptainer run -C -W${SLURM_TMPDIR} -B${WORKDIR}:/home -B /scratch -B /localscratch /scratch/zkileeg/Ecotype_Data/AGAT/agat_1.2.0.sif \
#	agat_sp_fix_overlaping_genes.pl --gff ${gff_in} -o ${gff_merged}  #merge overlapping loci
#apptainer run -C -W${SLURM_TMPDIR} -B${WORKDIR}:/home -B /scratch -B /localscratch /scratch/zkileeg/Ecotype_Data/AGAT/agat_1.2.0.sif \
#	agat_sp_keep_longest_isoform.pl -gff ${gff_merged} -o ${gff_rep} #find representative sequence
#apptainer run -C -W${SLURM_TMPDIR} -B${WORKDIR}:/home -B /scratch -B /localscratch /scratch/zkileeg/Ecotype_Data/AGAT/agat_1.2.0.sif \
#	agat_sp_extract_sequences.pl -f ${genome_in} -g ${gff_rep} -t cds -o ${cds_out} #output cds
#apptainer run -C -W${SLURM_TMPDIR} -B${WORKDIR}:/home -B /scratch -B /localscratch /scratch/zkileeg/Ecotype_Data/AGAT/agat_1.2.0.sif \
#	agat_sp_extract_sequences.pl -f ${genome_in} -g ${gff_rep} -p -o ${prot_out}   #output peptides

