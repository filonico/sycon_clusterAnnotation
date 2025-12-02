from samap.analysis import (get_mapping_scores)

from argparse import RawTextHelpFormatter
import pandas as pd
import sys, os, subprocess, pickle, argparse


subprocess.run("source ~/miniforge3/bin/activate SAMap_env", shell = True)


#####################
#     ARGUMENTS     #
#####################

parser = argparse.ArgumentParser(description = "Extract the SAMap mapping table from a pickle file.")
    
parser.add_argument("-p", "--input_pickle", required = True,
                    help = "Path to the pickle file to process.")

parser.add_argument("-o", "--output_dir", required = True,
                    help = "Path to the output directory.")

parser.add_argument("-n", "--n_top", default = "0",
                    help = "n_top parameter to use in the get_mapping_scores function. [default = 0]")

# check if the user gave no arguments, and if so then print the help
parser.parse_args(args = None if sys.argv[1:] else ["--help"])

args = parser.parse_args()


#######################
#     DEFINE I/Os     #
#######################

# define input file and output directory
output_dir = args.output_dir
input_pkl = args.input_pickle

# define other parameters which will be useful later on
ID = "originalVSremapped"

# create output dir if does not exist
if not os.path.isdir(output_dir):
    subprocess.run(f"mkdir -p {output_dir}",
                   shell = True)
    

##############################
#     GET MAPPING SCORES     #
##############################

# load the samap obj stored in the pickle file
with open(input_pkl, "rb") as file:
    samap_obj = pickle.load(file)

# create a dictionary for the cell cluster column name
cellCluster_dict = {"slacOriginal": "cell_type",
                    "slacRemapped": "seurat_clusters"}

# set the mapping score matrix name
matrix_filename = os.path.join(output_dir, ID + "_" +
                               args.n_top + "topCells_samapMappingTable.tsv")
scoringAln_filename = os.path.join(output_dir, ID + "_" +
                                   args.n_top + "topCells_samapScoringAln.tsv")
# get the mapping scores
D,MappingTable = get_mapping_scores(samap_obj, cellCluster_dict, n_top = int(args.n_top))

# save to files
MappingTable.to_csv(matrix_filename, index = True, sep = '\t')
D.to_csv(scoringAln_filename, index = False, sep = '\t')