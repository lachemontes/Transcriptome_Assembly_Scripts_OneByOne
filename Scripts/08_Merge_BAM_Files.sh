#!/bin/bash
#SBATCH -A naiss2023-5-461
#SBATCH -p core -n 16
#SBATCH -t 7-00:00:00
#SBATCH -J MergeBAM
#SBATCH --mail-type=All
#SBATCH --mail-user=zaide.montes_ortiz@biol.lu.se



# Load necessary modules or set necessary environment variables here if needed

module load bioinfo-tools
module load samtools/1.9



# Navigate to the directory containing sorted BAM files
cd /path/to/your/Bam_files/

# List all the sorted BAM files
bam_files=$(ls *.sorted.bam)

# Merge the BAM files
samtools merge merged_output.bam $bam_files

# Index the merged BAM file
samtools index merged_output.bam