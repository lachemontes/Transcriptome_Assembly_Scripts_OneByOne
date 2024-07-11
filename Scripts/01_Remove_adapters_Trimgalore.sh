#!/bin/bash
#SBATCH -A naiss2023-5-461
#SBATCH -p node
#SBATCH -t 2-00:00:00
#SBATCH -J Trimgalore
#SBATCH --mail-type=All
#SBATCH --mail-user=zaide.montes_ortiz@biol.lu.se
#SBATCH --array=1-96
#SBATCH --input=List_fastq_files.txt
#SBATCH --output=logs/trimgalore_%A_%a.out
#SBATCH --error=logs/trimgalore_%A_%a.err

# Name: Trimgalore_paired.sh
# Description: This script takes as input raw RNAseq data in fastq.gz listed in a txt file to remove adapters from Illumina sequencing.
# Author: Zaide Montes
# Institution: Lund University, Pheromone group to run in Rackham, uppmax
# Contact email: zk.montes10@gmail.com
# Date: Implemented on Mar 23, 2023

# Load necessary modules
# Do not change the cutadapt version used to avoid errors
module load bioinfo-tools
module load TrimGalore/0.6.1
module load FastQC/0.11.9
module load cutadapt/2.1

# Define input and output directories
fastqc_output_dir="/path/to/your/output/fastqc"
trim_galore_output_dir="path/to/your/output/Trimgalore"
log_dir="/path/to/your/output/Trimgalore_logs"

# Create necessary directories if they don't exist
mkdir -p $trim_galore_output_dir
mkdir -p $log_dir

# Get the sample name for this task from List.txt
sample=$(sed -n "${SLURM_ARRAY_TASK_ID} p" List_fastq_files.txt)

# Define input files
r1="${fastqc_output_dir}/${sample}_R1_001.fastq.gz"
r2="${fastqc_output_dir}/${sample}_R2_001.fastq.gz"

# Run TrimGalore with paired-end data, specifying to use illumina adapters, and run FastQC on the results
trim_galore -j 4 --illumina --paired "${r1}" "${r2}" --fastqc -o "${trim_galore_output_dir}"

# End of script
