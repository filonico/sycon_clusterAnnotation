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


#####################
#     FUNCTIONS     #
#####################

def add_speciesID_to_cluster_names(species_id):
    obs = samap_obj.samap.adata.obs.copy()
    
    # define column names
    cluster_col = f"{species_id}_N1.Int.Clusters"   # for sycon
    cell_type_col = f"{species_id}_cell_type"       # for any other species
    
    # define the criterion to match only corresponding lines
    mask = obs['species'] == species_id
    
    # add prefixes
    if cluster_col in obs.columns:
        obs.loc[mask, cluster_col] = (
            species_id + '_' + obs.loc[mask, cluster_col].astype(str)
        )
    
    if cell_type_col in obs.columns:
        obs.loc[mask, cell_type_col] = (
            species_id + '_' + obs.loc[mask, cell_type_col].astype(str)
        )
    
    # update the original dataframe
    samap_obj.samap.adata.obs = obs

#######################
#     DEFINE I/Os     #
#######################

# define input file and output directory
output_dir = args.output_dir
input_pkl = args.input_pickle
cell_cluster = args.cell_cluster
threshold = args.threshold

# define other parameters which will be useful later on
ID = os.path.basename(input_pkl).split("_")[0]              # ID of the SAMap run
species = [ID[:4], ID[-4:]]                                 # included species in the analysis
output_suffix = os.path.basename(input_pkl).split("_")[1]   # type of SAMap clustering

# load the table with species metadata
species_metadata = pd.read_table("00_input/species_metadata.tsv")

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
cellCluster_dict = species_metadata.set_index("speciesID")["cellCluster_annotation_name"].to_dict()

# add prefixes for the two species
for sp in species:
    add_speciesID_to_cluster_names(sp)

# find cluster specific markes
gpf = GenePairFinder(samap_obj, keys = cellCluster_dict)

# find gene pairs
if cell_cluster == "all":
    gene_pairs = gpf.find_all(align_thr = threshold)
else:
    gene_pairs = gpf.find_all(n = cell_cluster, align_thr = threshold)

# set the gene pair matrix output file name
matrix_gene_pairs = os.path.join(output_dir, ID + "_" +
                                 output_suffix + "_" +
                                 cell_cluster.replace("_", "") + "_samapGenePairs.tsv")

# save the gene pair matrix to a file
gene_pairs.to_csv(matrix_gene_pairs, index = False, sep = '\t')