#!/bin/env python3

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse, sys


##########################################
#     Define arguments of the script     #
##########################################

# Initialise the parser class
parser = argparse.ArgumentParser(description = "Get longest isoforms from fasta based on sequence identifiers.")

# Define some options/arguments/parameters
parser.add_argument("-i", "--input",
                    required = True,
                    help = "Input fasta file.")

parser.add_argument("-o", "--output",
                    help = "Output fasta file names. [default = \"inputName_noIso\"]")

# This line checks if the user gave no arguments, and if so then print the help
parser.parse_args(args = None if sys.argv[1:] else ["--help"])

# Collect the inputted arguments into a dictionary
args = parser.parse_args()

input_file = args.input
output_file = args.output


#############################
#     PROCESS THE FASTA     #
#############################

sequences_to_keep = {}
sequences_seen = []

for seq in SeqIO.parse(open(input_file), "fasta"):
    ID = seq.id
    SEQ = seq.seq
    if ID not in sequences_seen:
        sequences_to_keep[ID] = SEQ
        sequences_seen.append(ID)
    else:
        if len(SEQ) > len(sequences_to_keep[ID]):
            sequences_to_keep.update({ID: SEQ})

sequences_to_write = [SeqRecord(seq, id = id, description = "") for id, seq in sequences_to_keep.items()]        

with open(output_file, 'w') as handle:
    SeqIO.write(sequences_to_write, handle, 'fasta')
