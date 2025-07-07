#!/bin/bash

##############################
#     TIDY UP INPUT DATA     #
##############################

# tidy up raw data
# mind that data should be tar unzipped before
bash scripts/01_tidyup_data.sh


#####################
#     RUN BUSCO     #
#####################

mkdir 02_busco

# run busco on tidied proteomes
# REQUIRES: conda_envs/busco_env.yml
srun --partition=devel --job-name=physco --cpus-per-task=15 --time=02:00:00 --mem=4g --export=NONE --account=evassvis scripts/03_run_busco.sh

# summarise busco results into a tsv
grep "C:" 02_busco/*/short_summary*txt | sed -E 's/02_busco\///; s/_.+:\t/\t/; s/\t$//' > 02_busco/busco_result_summary.tsv

# compute species tree
# REQUIRES: conda_envs/phylo_env.yml
sbatch --job-name=physco --cpus-per-task=15 --time=04:00:00 --mem=4g --export=NONE --account=evassvis --mail-type=BEGIN,END,FAIL --mail-user=fn76@le.ac.uk scripts/04_compute_species_tree.sh

# gzip busco output directories
for i in 02_busco/*busco; do OUTNAME="$i".tar.gz; tar -cvf - "$i" | gzip -v -9 - > $OUTNAME && rm -rf $i; done


################################
#     RUN PAIRWISE DIAMOND     #
################################

mkdir 03_pairwise_diamond

# run all-vs-all diamond (exclude self comparison)
# REQUIRES: conda_envs/diamond_env.yml
srun --job-name=diamond --cpus-per-task=8 --time=04:00:00 --mem=0 --account=evassvis scripts/05_run_pairwise_diamond.sh


#############################
#     PREPARE H5AD FILES    #
############################

mkdir 04_preprocessed_scRNAseqs

# prepare .h5ad files
# REQUIRES: conda_envs/RSeurat_env.yml
Rscript scripts/06_prepare_h5ad_files.R


#####################
#     RUN SAMap     #
#####################

# run SAMap for sponges vs sponges
bash scripts/08_run_SAMap.sh Scil,Slac,Aque both 05b_SAMap_recodedSyconClusters | tee -a 05b_SAMap_recodedSyconClusters/Porifera_samap_both_leiden.log

# run SAMap for sycon vs placozoa
bash scripts/08_run_SAMap.sh Scil,Hhon,HH23,Tadh,TrH2 both 05b_SAMap_recodedSyconClusters | tee -a 05b_SAMap_recodedSyconClusters/ScilPlacozoa_samap_both_leiden.log

# run SAMap for spongilla vs placozoa
bash scripts/08_run_SAMap.sh Slac,Hhon,HH23,Tadh,TrH2 both 05b_SAMap_recodedSyconClusters | tee -a 05b_SAMap_recodedSyconClusters/SlacPlacozoa_samap_both_leiden.log

# run SAMap for sycon vs mnemiopsis
bash scripts/08_run_SAMap.sh Scil,Mlei pairwise 05b_SAMap_recodedSyconClusters | tee -a 05b_SAMap_recodedSyconClusters/ScilMlei_samap_pairwise_leiden.log

# run SAMap for sycon vs cnidaria
bash scripts/08_run_SAMap.sh Scil,Nvec,Spis,Hvul,Xesp both 05b_SAMap_recodedSyconClusters | tee -a 05b_SAMap_recodedSyconClusters/ScilCnidaria_samap_both_leiden.log


################################
#     Get SAMap statistics     #
################################

# get mapping scores
for i in 05b_SAMap_recodedSyconClusters/*pkl; do python scripts/09_get_SAMap_mappingTables.py -p $i -o 05b_SAMap_recodedSyconClusters/01_mapping_scores -n 0; done

# get gene pairs for Placozoa and Sycon cluster 14
for i in 05b_SAMap_recodedSyconClusters/{HhonScil,HH23Scil,ScilTadh,ScilTrH2}*pkl; do python scripts/10_get_SAMap_genePairs.py -p $i -o 05b_SAMap_recodedSyconClusters/02_gene_pairs -c Scil_14 -t 0.2; done

# get gene pairs for Placozoa and all Sycon clusters
for i in 05b_SAMap_recodedSyconClusters/{HhonScil,HH23Scil,ScilTadh,ScilTrH2}*pkl; do python scripts/10_get_SAMap_genePairs.py -p $i -o 05b_SAMap_recodedSyconClusters/02_gene_pairs -t 0.2; done
