#!/bin/bash
#SBATCH -A snic2022-5-454
#SBATCH -p core -n 12
#SBATCH -t 3-00:00:00
#SBATCH -J Transdecoder
#BATCH --mail-type=all
#SBATCH --mail-user=zaide.montes_ortiz@biol.lu.se


module load bioinfo-tools
module load TransDecoder/5.7.0


TransDecoder.LongOrfs -t /path/to/your/transcripts.fasta
TransDecoder.Predict -t /path/to/your/transcripts.fasta