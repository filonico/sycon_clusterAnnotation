#!/bin/env python3


# Given a list file of SRA accession numbers, this script is intented to download the corresponding reads from NCBI.
# REQUIRED SOFTWARES: sra-toolkit, fastqc
#
# The final output would be a directory structured as follow:
#
# ./
# └── your_output_dir/
#     ├── 01_fastq/
#     |   └── {results of fastqc analysis}
#     ├── SRRXXXXXX1/
#     |   ├── SRRXXXXXX1_1.fastq.gz
#     |   └── SRRXXXXXX1_2.fastq.gz
#     ├── SRRXXXXXX2/
#     |   ├── SRRXXXXXX2_1.fastq.gz
#     |   └── SRRXXXXXX2_2.fastq.gz
#     ...
#     └── SRRXXXXXXN/
#         ├── SRRXXXXXXN_1.fastq.gz
#         └── SRRXXXXXXN_2.fastq.gz
#
#
# Written by:   Filippo Nicolini
# Last updated: 18/10/2023


import subprocess, argparse, sys, os


##########################################
#     Define arguments of the script     #
##########################################

# Initialise the parser class
parser = argparse.ArgumentParser(description = "Download reads from NCBI through the sra-tool and perform quality check.")

# Define some options/arguments/parameters
parser.add_argument("-i", "--input",
                    required = True,
                    help = "A list of SRA accession numbers.")

parser.add_argument("-o", "--output_dir",
                    help = "Name of the output directory.",
                    default = "01_raw_reads")

# This line checks if the user gave no arguments, and if so then print the help
parser.parse_args(args = None if sys.argv[1:] else ["--help"])

# Collect the inputted arguments into a dictionary
args = parser.parse_args()


output_dir = args.output_dir


############################
#     Define functions     #
############################

# Function to download the SRA file, given the SRA accession number
def download_sra(sra):
    try:
        # Download SRA
        prefetch_process = subprocess.run(f"prefetch -O {output_dir}/ {sra}",
                                          shell = True,
                                          capture_output = True,
                                          text = True)
        
    except subprocess.CalledProcessError as err:
        print("An error occured:", err.stderr)


# Function to download the fastq file, given the SRA file
def download_fastq(sra_file):
    try:
        # Download fastqs
        fastqdump_process = subprocess.run(f"fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files {sra_file} -O {fastq_output_dir}",
                                           shell = True,
                                           capture_output = True,
                                           text = True)
    
    except subprocess.CalledProcessError as err:
        print("An error occured:", err.stderr)
        

# Function to perform the quality check, given the fastq file
def quality_check(fastq_file):
    try:
        fastqc_process = subprocess.run(f"fastqc {fastq_file} -o {output_dir}/01_fastqc -f fastq",
                                        shell = True,
                                        capture_output = True,
                                        text = True)

    except subprocess.CalledProcessError as err:
        print("An error occured:", err.stderr)



#------------------------------------------------------------------------------------------


#########################
#     Read SRA list     #
#########################

if not os.path.isdir(args.output_dir):
    print()
    print(f"Creating output directory in {args.output_dir}/")
    subprocess.run(f"mkdir {args.output_dir}", shell = True)

# Read in SRA list and store into a list object
SRA_list = []
with open(args.input) as input_SRA:
    for line in input_SRA.readlines():
        SRA_list.append(line.strip())

print()
print(f"Read {args.input}: {len(SRA_list)} accession numbers found")
print()


####################################################
#     Download reads and perform quality check     #
####################################################

for SRA in SRA_list:

    # Download SRA
    print(f"-- {SRA} --")
    print("  Retrieving SRA files...")
    download_sra(SRA)

    # Download fastq
    SRA_file = output_dir + "/" + SRA + "/*sra*"
    fastq_output_dir = output_dir + "/" + SRA
    print("  Retrieving fastq files...")
    download_fastq(SRA_file)

    # Quality check
    FASTQ_file = fastq_output_dir + "/*fastq"
    print("  Checking read quality...")
    quality_check(FASTQ_file)

    # Gzip fastq files
    print("  Gzipping fastq files...")
    subprocess.run(f"gzip -9 {fastq_output_dir}/*fastq",
                    shell = True)
    
    # Remove sra files
    print("  Removing the sra file...")
    subprocess.run(f"rm {SRA_file}",
                    shell = True)
    
    print("Done")
    print()
