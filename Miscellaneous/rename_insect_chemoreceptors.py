"""
Script to rename FASTA files of insect chemoreceptors from RefSeq.
The script reads FASTA files containing sequences and their descriptions,
and renames each sequence according to the specified format:
accessionnumber_Abbreviation_GR_Number.

Requirements
Make sure to have the Biopython library installed:
pip install biopython


Usage:
    python rename_insect_chemoreceptors.py <file1.fasta> [<file2.fasta> ...]

Requirements:
    - The script requires the Biopython library to be installed.
      Install it using: pip install biopython
"""

from Bio import SeqIO
import re  
import sys
import os

# Ensure that at least one file is provided as an argument
if len(sys.argv) < 2:
    print("Usage: python rename_insect_chemoreceptors.py <file1.fasta> [<file2.fasta> ...]")
    sys.exit(1)

# Process each input file provided through command line arguments
for fasta_file in sys.argv[1:]:
    filename = os.path.splitext(fasta_file)[0]  # Get the filename without extension
    output_file = "Fix_" + filename + ".fasta"  # Define output filename

    with open(fasta_file, "r") as f_open, open(output_file, "a+") as specie_file:
        for rec in SeqIO.parse(f_open, "fasta"):
            description = rec.description
            
            # Extract the number from the "GR" (Gustatory Receptor)
            rex1 = 'GR(\d+)'  # Capture the number that follows 'GR'
            protein_number_matches = re.search(rex1, description)
            protein_number = protein_number_matches.group(1) if protein_number_matches else ""

            # Extract the number that follows "gustatory receptor"
            gustatory_receptor_match = re.search(r'gustatory receptor (\d+)', description)
            if gustatory_receptor_match:
                protein_number += gustatory_receptor_match.group(1)  # Append the found number

            # Extract the species name
            specie = description[description.find("[") + 1:description.find("]")]
            short_specie = specie.split()[0][0] + specie.split()[1][0:3]  # Example: Aedes aegypti -> Aaeg
            
            # Construct the new ID
            new_id = f"{rec.id}_{short_specie}_GR{protein_number}"  # Desired format: accession_abbreviation_GR_number
            
            # Write to the output file
            specie_file.write(f">{new_id}\n{rec.seq}\n")
