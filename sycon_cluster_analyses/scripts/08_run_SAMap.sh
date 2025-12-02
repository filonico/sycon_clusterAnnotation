#!/bin/bash

# $1 = comma-separeted list of species to analyse
# $2 = type of SAMap analysis (pairwise, stitched, or both)
# $3 = output directory

source ~/miniforge3/bin/activate SAMap_env

python scripts/07_run_SAMap.py -s $1 -a $2 -o $3
