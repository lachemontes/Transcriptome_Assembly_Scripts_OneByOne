#!/bin/bash
#SBATCH -A naiss2023-5-461
#SBATCH -p core -n 16
#SBATCH -t 7-00:00:00
#SBATCH -J HISAT_Ticks_Index
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zaide.montes_ortiz@biol.lu.se
#SBATCH --output=logs/hisat_%A_%a.out
#SBATCH --error=logs/hisat_%A_%a.err

# Name: HISAT2_Index.sh
# Description: This script builds a HISAT2 index for the reference genome.
# Author: Zaide Montes
# Institution: Lund University, Pheromone group to run in Rackham, uppmax
# Contact email: zk.montes10@gmail.com
# Date: Implemented on Nov 7, 2022

# Load necessary modules
module load bioinfo-tools HISAT2/2.2.1

# Define directories
genome_file="/path/to/your/genome_file.fna"
index_output_dir="/path/to/your/HISAT/Genome_Index"
log_dir="path/to/your/HISAT_logs"

# Create necessary directories if they don't exist
mkdir -p $index_output_dir
mkdir -p $log_dir

# Run HISAT2 build to create genome index
hisat2-build -p 16 "$genome_file" "$index_output_dir/Genome_Index"

# End of script
