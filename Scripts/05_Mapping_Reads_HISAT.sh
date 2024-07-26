#!/bin/bash
#SBATCH -A snic2022-5-454
#SBATCH -p node
#SBATCH -t 12:00:00
#SBATCH -J Fem_HISAT_Mapping
#SBATCH --array=1-96
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zaide.montes_ortiz@biol.lu.se
#SBATCH --output=logs/hisat_%A_%a.out
#SBATCH --error=logs/hisat_%A_%a.err

# Name: HISAT2_Mapping.sh
# Description: This script maps RNAseq reads, after adapter trimming with TrimGalore, to a reference genome using HISAT2.
# Author: Zaide Montes
# Institution: Lund University, Pheromone group to run in Rackham, uppmax
# Contact email: zk.montes10@gmail.com
# Date: Implemented on Nov 7, 2022

# Load necessary modules
module load bioinfo-tools
module load HISAT2/2.2.1

# Define input and output directories
trim_galore_output_dir="path/to/your/output/Trimgalore"
hisat_output_dir="path/to/your/output/HISAT"
log_dir="/path/to/your/output/HISAT_logs"
genome_index="/path/to/your/output/HISAT/Genome_Index"

# Create necessary directories if they don't exist
mkdir -p $hisat_output_dir
mkdir -p $log_dir

# Get the sample name for this task from List.txt
sample=$(sed -n "${SLURM_ARRAY_TASK_ID} p" List.txt)

# Define input files
r1="${trim_galore_output_dir}/${sample}_R1_001_val_1.fq.gz"
r2="${trim_galore_output_dir}/${sample}_R2_001_val_2.fq.gz"

# Define output file
output="${hisat_output_dir}/${sample}.sam"

# Run HISAT2
hisat2 --dta -p 8 -x "${genome_index}" -1 "${r1}" -2 "${r2}" -S "${output}"

# Merge HISAT *.err files into a consolidated log file
tail -n +1 logs/*.err > logs/concatenated_hisat_errors.txt

# End of script
