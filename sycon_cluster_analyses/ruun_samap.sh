#!/bin/bash

source ~/miniforge3/bin/activate samap_env

python scripts/07NEW_run_SAMap.py -s Aque,Scil,Slac,Hvul,Nvec,Spis,Xesp -a both -i 04_preprocessed_scRNAseqs/ -d 03_pairwise_diamond/ -n 6 -o 15NEW_SAMap_cnidaria/
