#!/bin/bash

while getopts n:s:y: flag
do
    case "${flag}" in
        n) NCORES=${OPTARG};;
        s) STEP=${OPTARG};;
        y) YAML_PATH=${OPTARG};;
    esac
done

PATH_TO_CODE="/lustre/alice3/data/evassvis/fn76/SOFTWARES/Bonsai-data-representation"

srun python -m mpi4py ${PATH_TO_CODE}/bonsai/bonsai_main.py --config_filepath ${YAML_PATH} --step ${STEP}
