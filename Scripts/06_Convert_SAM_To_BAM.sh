#!/bin/bash
#SBATCH -A naiss2023-5-461
#SBATCH -p node
#SBATCH -t 7-00:00:00
#SBATCH -J SAMToBAM
#SBATCH --mail-type=All
#SBATCH --mail-user=zaide.montes_ortiz@biol.lu.se
#SBATCH --array=1-96

# Load necessary modules or set necessary environment variables here if needed

module load bioinfo-tools
module load samtools/1.9

filename=$(sed -n "${SLURM_ARRAY_TASK_ID}p" List.txt)
input_file="/path/to/your/${filename}.sam"
output_file="/path/to/your/output/Bam/${filename}.bam"
sorted_output_file="/path/to/your/output/Bam/${filename}.sorted.bam"

# Convert SAM to BAM
samtools view -bS $input_file > $output_file

# Sort the BAM file
samtools sort $output_file -o $sorted_output_file

# Index the sorted BAM file if needed
samtools index $sorted_output_file