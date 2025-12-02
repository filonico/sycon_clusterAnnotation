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

mapping_dir = "./04_SAMap/01_pairwise_diamond/"
h5ads_dir = "./04_SAMap/02_prepared_h5ads/"
output_dir = "./04_SAMap/"


########################
#     IMPORT FILES     #
########################

# define the dictionary  where spIDs and path to h5ad files will be stored
data = {}

# populate the dictionary
for filepath, dirnames, filenames in os.walk(h5ads_dir):
    for filename in filenames:
        print("Loading raw data...")
        # extract the species ID
        spID = filename[0:12]
        # get the path to the file
        h5ad_file = os.path.join(filepath, filename)
        # add the dict entry
        data[spID] = h5ad_file

print("\n")

print(data)


###################################
#     RUN THE SAMAP ALGORITHM     #
###################################

resolution_dict = {key: 3 for key in data}
speciesTrue_dict = {key: True for key in data}

sam_mapping = SAMAP(data,
                    f_maps = mapping_dir,
                    resolutions = resolution_dict,
                    save_processed = False)

sam_mapping.run(pairwise = True,
                neigh_from_keys = speciesTrue_dict,
                NUMITERS = 6)

with open(os.path.join(output_dir, "slac_integrated_originalVSremapped.pkl"), "wb") as file:
    pickle.dump(sam_mapping, file)

sam_mapping.samap.save_anndata(os.path.join(output_dir, "slac_integrated_originalVSremapped.h5ad"), "wb")