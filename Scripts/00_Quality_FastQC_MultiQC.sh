#!/bin/bash
#SBATCH -A snic2022-5-454
#SBATCH -p core -n 8
#SBATCH -t 12:00:00
#SBATCH -J Fastqc
#SBATCH --mail-type=All
#SBATCH --mail-user=zaide.montes_ortiz@biol.lu.se
#SBATCH --array=1-96
#SBATCH --input=List_fastq_files.txt
#SBATCH --output=logs/fastqc_%A_%a.out
#SBATCH --error=logs/fastqc_%A_%a.err

# Directories for input and output files
INPUT_DIR="/path/to/your/input/files"
OUTPUT_DIR="/path/to/your/output/fastqc"
MULTIQC_DIR="/path/to/your/output/multiqc"
LOG_DIR="/path/to/your/output/logs"

# Create necessary directories if they don't exist
mkdir -p $OUTPUT_DIR
mkdir -p $MULTIQC_DIR
mkdir -p $LOG_DIR

# Load FastQC module
module load bioinfo-tools
module load FastQC/0.11.9
module load MultiQC/1.12

# Get the input file for this task by reading the corresponding line in List_fastq_files.txt
input_file=$(sed -n "$SLURM_ARRAY_TASK_ID p" List_fastq_files.txt)

# Run FastQC on the input file, output results to the designated directory
fastqc -o $OUTPUT_DIR $INPUT_DIR/$input_file

# Run MultiQC and output the summary to the designated directory
multiqc -o $MULTIQC_DIR $OUTPUT_DIR

# End of script
