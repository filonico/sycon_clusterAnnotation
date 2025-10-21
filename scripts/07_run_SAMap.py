#!/usr/bin/env python

import subprocess

subprocess.run("source ~/miniforge3/bin/activate SAMap_env", shell = True)

from samap.mapping import SAMAP
from samalg import SAM
import pandas as pd
import sys, os, pickle, argparse
from argparse import RawTextHelpFormatter


##################
#     INPUTS     #
##################

mapping_dir = './03_pairwise_diamond/'

# load the table with species metadata
species_metadata = pd.read_table("00_input/species_metadata.tsv")


#####################
#     FUNCTIONS     #
#####################

# check if input species IDs are valid or not (based on species_metadata file)
def validate_species(species_list):
    species = sorted([s.strip() for s in species_list.split(",")])
    if not set(species).issubset(species_metadata["speciesID"]):
        invalid_species = set(species).difference(species_metadata["speciesID"])
        raise argparse.ArgumentTypeError(f"Invalid species ID(s): {', '.join(invalid_species)}")
    if len(species) <= 1:
        raise argparse.ArgumentTypeError("You must select at least two species.")
    return species

# create samap object with pre-computed cell clusters
def create_samap_with_precomputed_clusters(input):
    # create samap object
    samap_object = SAMAP(input,
                        f_maps = mapping_dir,
                        keys = cellCluster_dict,
                        save_processed = False)
    return(samap_object)

# create samap object with samap-computed cell clusters
def create_samap_with_leiden_clustering(input):
    # assign the default resolution value to all species
    resolution_dict = {key: 3 for key in species_to_process}
    # create samap object
    samap_object = SAMAP(input,
                         f_maps = mapping_dir,
                         resolutions = resolution_dict,
                         save_processed = False)
    return(samap_object)

# update .obs in the anndata object with full species names
def update_sam_mapping_obs(samap_object):
    samap_object.samap.adata.obs = samap_object.samap.adata.obs.merge(species_metadata[["speciesID", "species_full"]],
                                                                      left_on = "species",
                                                                      right_on = "speciesID",
                                                                      how = "inner").drop(columns = "speciesID")
    
def save_processed_sam_objects(samap_objects, file_annotation):
    for spID in samap_objects.sams:
        samap_objects.sams[spID].save_anndata(os.path.join(input_dir, spID + '_cellFiltered_' + file_annotation + '.h5ad'))


#####################
#     ARGUMENTS     #
#####################

parser = argparse.ArgumentParser(description = "Run SAMap on pre-filtered single-cell h5ad files.")
    
parser.add_argument("-s", "--species",
                    type = validate_species, required = True,
                    help = "Comma-separated list of species IDs (e.g., species1,species2,species3)")

parser.add_argument("-a", "--analysis",
                    choices = ["pairwise", "stitched", "both"], default = "both",
                    help = "Type of analysis. Choose among \"pairwise\", \"stitched\", and \"both\". [default: both]")

parser.add_argument("-d", "--data",
                    choices = ["raw", "preprocessed"], default = "raw",
                    help = "Type of data to load. Choose between \"raw\" (for raw counts) and \"preprocessed\" (for preprocessed SAM files). Must be in h5ad files. [default: raw]")

parser.add_argument("-c", "--clustering",
                    choices = ["precomputed", "leiden"], default = "leiden",
                    help = "Clustering method. Choose between \"precomputed\" (if you already know the clusters) and \"leiden\" (to allow SAMap do the leiden clustering). [default : leiden]")

parser.add_argument("-i", "--input_dir",
                    default = "./04_preprocessed_scRNAseqs/",
                    help = "Input directory where h5ad files are stored. [default : ./04_preprocessed_scRNAseqs/]")


parser.add_argument("-o", "--output_dir", required = True,
                    help = "Output directory.")


# check if the user gave no arguments, and if so then print the help
parser.parse_args(args = None if sys.argv[1:] else ["--help"])

args = parser.parse_args()

type_of_analysis_dict = {
    "pairwise": "SAMap true pairwise",
    "stitched": "SAMap stitched species",
    "both": "SAMap pairwise + stitched"
}
input_data_type_dict = {
    "raw": "raw counts matrices",
    "preprocessed": "pre-processed SAM matrices"
}
clustering_method_dict = {
    "precomputed": "pre-computed cell clusters",
    "leiden": "SAMap clusters with resolution 3"
}


##########################
#     ANALYSIS RECAP     #
##########################

input_dir = args.input_dir
output_dir = args.output_dir

# create output dir if does not exist
if not os.path.isdir(output_dir):
    subprocess.run(f"mkdir -p {output_dir}",
                   shell = True)

print("\n#####   ANALYSIS RECAP   #####")
print(f"Species to process: {', '.join(args.species)}")
print(f"File to load: {input_data_type_dict[args.data]}")
print(f"Type of analysis: {type_of_analysis_dict[args.analysis]}")
print(f"Cell-clustering method: {clustering_method_dict[args.clustering]}")
print(f"Input directory: {input_dir}")
print(f"Output directory: {output_dir}")
print("##############################\n")

sys.stdout.flush()

###############################
#     SET UP THE ANALYSIS     #
###############################

species_to_process = args.species

if args.clustering == "precomputed":
    output_suffix = 'preComputedClusters'
elif args.clustering == "leiden":
    output_suffix = 'leiden3Clusters'

# sys.exit("\nNothing more to do for now :-) we are tetsing the code. Exiting...\n")

# subset the dataset to only species that needs to be analysed
species_metadata = species_metadata[species_metadata["speciesID"].isin(species_to_process)]


########################
#     IMPORT FILES     #
########################

# define the dictionary  where spIDs and path to h5ad files will be stored
data = {}

# populate the dictionary
for filepath, dirnames, filenames in os.walk(input_dir):
    for filename in filenames:
        # load only files of chosen species
        if any(species in filename for species in species_to_process):
            # if raw data have been chosen, then load them
            if args.data == "raw" and "Clusters" not in filename:
                print("Loading raw data...")
                # extract the species ID
                spID = filename[0:4]
                # get the path to the file
                h5ad_file = os.path.join(filepath, filename)
                # add the dict entry
                data[spID] = h5ad_file
            # if preprocessed SAM data have been chosen, then load them
            elif args.data == "precomputed" and "Clusters" in filename:
                print("Loading pre-processed data...")
                # extract the species ID
                spID = filename[0:4]
                # create an empty SAM object
                preprocessed_file = SAM()
                # load the preprocessed file as a SAM object
                preprocessed_file.load_data(os.path.join(filepath, filename.split('_')[0] + '_cellFiltered_' + output_suffix + '.h5ad'))
                # add the dict entry
                data[spID] = preprocessed_file

print("\n")

print(data)

# create a dictionary for the name of cell cluster annotations
cellCluster_dict = species_metadata.set_index("speciesID")["cellCluster_annotation_name"].to_dict()

sys.stdout.flush()

###################################
#     RUN THE SAMAP ALGORITHM     #
###################################

# create a dictionary with True booleans to be used in the neigh_from_keys SAMap algorithm
speciesTrue_dict = {key: True for key in species_metadata["speciesID"]}

# STITCHED SPECIES ANALYSIS
if args.analysis == "stitched" or args.analysis == "both":

    # PRE-COMPUTED CELL CLUSTERS
    if args.clustering == "precomputed":
        # use the pre-computed cell clusters
        sam_mapping = create_samap_with_precomputed_clusters(data)
        save_processed_sam_objects(sam_mapping, output_suffix)
        # run the SAMAP algorithm
        sam_mapping.run(pairwise = True, neigh_from_keys = speciesTrue_dict,
                        NUMITERS = 6)

    # SAMAP CELL CLUSTERING
    elif args.clustering == "leiden":
        # use the SAMap clustering method, assuming a resolution parameter of three for all species
        sam_mapping = create_samap_with_leiden_clustering(data)
        save_processed_sam_objects(sam_mapping, output_suffix)
        # run the SAMAP algorithm
        sam_mapping.run(pairwise = True,
                        neigh_from_keys = speciesTrue_dict,
                        NUMITERS = 6)

    # update .obs in the anndata object with full species names
    update_sam_mapping_obs(sam_mapping)

    sys.stdout.flush()

    # save the computed SAMAP object
    with open(os.path.join(output_dir, "".join(species_to_process) + "_" + output_suffix + "_samap.pkl"), "wb") as file:
        pickle.dump(sam_mapping, file)

# TRUE PAIRWISE ANALYSIS
if args.analysis == "pairwise" or args.analysis == "both":
    analysis_pairs = [[x, y] for i, x in enumerate(species_to_process) for y in species_to_process[i + 1:]]

    for pair in analysis_pairs:
        # subset the data on the considered pair
        data_pairs = {species: data[species] for species in pair}

        # PRE-COMPUTED CELL CLUSTERS
        if args.clustering == "precomputed":
            # use the pre-computed cell clusters
            sam_mapping = create_samap_with_precomputed_clusters(data_pairs)
            save_processed_sam_objects(sam_mapping, output_suffix)
            # run the SAMAP algorithm
            sam_mapping.run(pairwise = True, neigh_from_keys = speciesTrue_dict,
                            NUMITERS = 6)

        # SAMAP CELL CLUSTERING
        elif args.clustering == "leiden":
            # use the SAMap clustering method, assuming a resolution parameter of three for all species
            sam_mapping = create_samap_with_leiden_clustering(data_pairs)
            save_processed_sam_objects(sam_mapping, output_suffix)
            # run the SAMAP algorithm
            sam_mapping.run(pairwise = True,
                            neigh_from_keys = speciesTrue_dict,
                            NUMITERS = 6)

        # update .obs in the anndata object with full species names
        update_sam_mapping_obs(sam_mapping)

        # save the computed SAMAP object
        with open(os.path.join(output_dir, "".join(pair) + "_" + output_suffix + "_samap.pkl"), "wb") as file:
            pickle.dump(sam_mapping, file)

        sys.stdout.flush()
