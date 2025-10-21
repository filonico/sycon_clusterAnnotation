#!/bin/bash

source ~/miniforge3/bin/activate bonsai_env

for i in 08_bonsai/S2*yaml; do FULLPATH="$(realpath $i)"; python3 ../SOFTWARES/Bonsai-data-representation/bonsai/bonsai_main.py --config_filepath $FULLPATH --step all; done
