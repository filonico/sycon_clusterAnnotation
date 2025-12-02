#!/bin/bash

source ~/miniforge3/bin/activate phySCO_env

# create directory to link proteomes
mkdir -p 02_busco/01_species_tree_aa/00_input

# create symlinks for proteomes
for i in 02_busco/*faa*; do LINK="$(realpath $i)" && ln -s $LINK 02_busco/01_species_tree_aa/00_input/; done

# compute specie tree through phySCO
python3 precompiled_softwares/phySCO/phySCO.py -i 02_busco/01_species_tree_aa/00_input/ -o 02_busco/01_species_tree_aa/01_ML_tree -g 100 -m
