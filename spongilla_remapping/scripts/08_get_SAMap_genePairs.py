from samap.analysis import (GenePairFinder)

from argparse import RawTextHelpFormatter
import pandas as pd
import sys, os, subprocess, pickle, argparse


subprocess.run("source ~/miniforge3/bin/activate SAMap_env", shell = True)


#####################
#     ARGUMENTS     #
#####################

parser = argparse.ArgumentParser(description = "Extract the gene pairs from SAMap results. NOTE that this script will work only with species pairs at the moment.")
    
parser.add_argument("-p", "--input_pickle", required = True,
                    help = "Path to the pickle file to process.")

parser.add_argument("-o", "--output_dir", required = True,
                    help = "Path to the output directory.")

parser.add_argument("-c", "--cell_cluster", default = "all",
                    help = "Name of the cell cluster of which to get gene pairs. Please prepend the species ID to the name (e.g., \"Scil_pinacocytes\"). [default = all]")

parser.add_argument("-t", "--threshold", default = 0.25, type = float, 
                    help = "Mapping score threshold. [default = 0.25]")


# check if the user gave no arguments, and if so then print the help
parser.parse_args(args = None if sys.argv[1:] else ["--help"])

args = parser.parse_args()


#######################
#     DEFINE I/Os     #
#######################

# define input file and output directory
output_dir = args.output_dir
input_pkl = args.input_pickle
cell_cluster = args.cell_cluster
threshold = args.threshold

# define other parameters which will be useful later on
ID = "originalVSremapped"                                   # ID of the SAMap run
#species = [ID[:4], ID[-4:]]                                 # included species in the analysis

# create output dir if does not exist
if not os.path.isdir(output_dir):
    subprocess.run(f"mkdir -p {output_dir}",
                   shell = True)
    

##########################
#     GET GENE PAIRS     #
##########################

# load the samap obj stored in the pickle file
with open(input_pkl, "rb") as file:
    samap_obj = pickle.load(file)

# create a dictionary for the cell cluster column name
cellCluster_dict = {"slacOriginal": "cell_type",
                    "slacRemapped": "seurat_clusters_2"}

# find cluster specific markes
gpf = GenePairFinder(samap_obj, keys = cellCluster_dict)

# find gene pairs
if cell_cluster == "all":
    gene_pairs = gpf.find_all(align_thr = threshold)
else:
    gene_pairs = gpf.find_all(n = cell_cluster, align_thr = threshold)

# set the gene pair matrix output file name
matrix_gene_pairs = os.path.join(output_dir, ID + "_" +
                                 cell_cluster.replace("_", "") + "_samapGenePairs.tsv")

# save the gene pair matrix to a file
gene_pairs.to_csv(matrix_gene_pairs, index = False, sep = '\t')
