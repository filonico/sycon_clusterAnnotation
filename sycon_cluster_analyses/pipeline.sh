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
srun --partition=devel --job-name=busco --cpus-per-task=15 --time=02:00:00 --mem=4g --export=NONE --account=evassvis scripts/03_run_busco.sh

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


####################################
#     GENE EXPRESSION ANALYSIS     #
####################################

# here goes the code to check for expression of specific genes/pathways


##################
#     BONSAI     #
##################

# run sanity on sycon UMI counts
# REQUIRES: conda_envs/sanity_env.yml
for i in 08_bonsai/00_umis_metadata/*_raw_UMImatrix.tsv; do SAMPLE=$(echo $i | sed -E 's/^.+\///; s/_raw.+$//'); OUTDIR=$(echo 08_bonsai/01_sanity_extOutput/"$SAMPLE"_sanity_extOutput); mkdir $OUTDIR && ../SOFTWARES/Sanity/bin/Sanity -f $i -n 15 -e true -max_v true -d $OUTDIR; done

# create bonsai config file (on spartacus)
for i in 08_bonsai/01_sanity_extOutput/*; do FULLPATH="$(realpath $i)"; NAME="$(basename $i | awk -F "_" '{print $1}')"; DIRPATH="$(realpath $(dirname $i | awk -F "/" '{print $1}'))"; python3 ../../SOFTWARES/Bonsai-data-representation/bonsai/create_config_file.py --new_yaml_path "$DIRPATH"/"$NAME"_bonsai_config.yaml --dataset $NAME --data_folder $FULLPATH --verbose True --results_folder $DIRPATH/02_"$NAME"_results --input_is_sanity_output True; done

# run bonsai (on spartacus)
for i in 08_bonsai/S2*yaml; do FULLPATH="$(realpath $i)"; python3 ../../SOFTWARES/Bonsai-data-representation/bonsai/bonsai_main.py --config_filepath $FULLPATH; done

# preprocess before visualising with bonsai-scout (on spartacus)
python3 ../../SOFTWARES/Bonsai-data-representation/bonsai_scout/bonsai_scout_preprocess.py --results_folder /DataDrives/Drive2/Filippo/ANALYSES/sycon_bonsai/08_bonsai/02_sycon_SCT_results/ --annotation_path /DataDrives/Drive2/Filippo/ANALYSES/sycon_bonsai/00_input/sycon_metadata.tsv --take_all_genes False --config_filepath ''

# run bonsai scout (on spartacus)
python3 ../../SOFTWARES/Bonsai-data-representation/bonsai_scout/run_bonsai_scout_app.py --results_folder /DataDrives/Drive2/Filippo/ANALYSES/sycon_bonsai/08_bonsai/02_sycon_SCT_results/ --settings_filename bonsai_vis_settings.json --port 1234


###############################
#     SEQUENCE ANNOTATION     #
###############################

bash /lustre/alice3/data/evassvis/fn76/SOFTWARES/InterProScan/interproscan-5.75-106.0/interproscan.sh -i 01_proteomes/Scil_ol.faa -goterms -b 09_gene_annotation/scil_proteome_interproscan


#################################
#     PERFORM GO ENRICHMENT     #
#################################

# GO terms were annotated from the Sycon proteome with the OMA Web Server

# get the list of cluster markers and the gene universe
Rscript scripts/17_get_markers_forGOenrich.R

# get GO annotation for each cluster marker and the gene universe
for i in 10_GO_enrichment/*ls; do grep -wf $i 09_gene_annotation/GOterms_OMA.tsv > "${i%%.*}"_GOterms.tsv; done

# perform GO enrichment for each cluster
for i in 10_GO_enrichment/cluster*tsv; do Rscript scripts/18_perform_GOenrich.R 10_GO_enrichment/geneUniverse_GOterms.tsv $i "${i%%.*}"_; done


###########################
#     KEGG ENRICHMENT     #
###########################

# selected organisms (eukaryotes + porifera + placo + some cnidarians): hsa, mmu, rno, dre, dme, cel, ath, sce, ago, cal, spo, ecu, pfa, cho, ehi, eco, nme, hpy, bsu, lla, mge, mtu, syn, aae, mja, ape, aqu, tad, nve, epa, adf, amil, pdam, spis, dgt, hmg
# assigned method: BBH

mkdir 11_KEGG_enrichment

# get the list of cluster markers
for i in 10_GO_enrichment/cluster*.ls; do ln -s $(realpath $i) 11_KEGG_enrichment/; done

# get gene universe annotation
grep -wf 10_GO_enrichment/geneUniverse.ls 09_gene_annotation/KOterms_kaas.tsv > 11_KEGG_enrichment/geneUniverse_KOterms.tsv

# perform KO enrichment for each cluster
for i in 11_KEGG_enrichment/cluster*ls; do Rscript scripts/19_perform_KEGGenrich.R 11_KEGG_enrichment/geneUniverse_KOterms.tsv $i "${i%%.*}"_KOenrich.tsv && echo done_$i; done


###################
#     hdWGCNA     #
###################

# CODE TO RUN hdWGCNA
Rscript scripts/20_hdWGCNA.R

# get GO annotation for each module
for i in 12_hdWGCNA/*ls; do grep -wf $i 09_gene_annotation/GOterms_OMA.tsv > "${i%%.*}"_GOterms.tsv; done

# perform GO enrichment for each module
for i in 12_hdWGCNA/*tsv; do Rscript scripts/18_perform_GOenrich.R 10_GO_enrichment/geneUniverse_GOterms.tsv $i "${i%%.*}"_; done

# perform KO enrichment for each module
for i in 12_hdWGCNA/*ls; do Rscript scripts/19_perform_KEGGenrich.R 11_KEGG_enrichment/geneUniverse_KOterms.tsv $i "${i%%.*}"_KOenrich.tsv && echo done_$i; done


######################################
#     RECLUSTER THE CENTRAL BLOB     #
######################################

mkdir -p 13_recluster_blob/{01_onlyBlob_originalClusters,02_onlyBlob_newClusters,03_hdWGCNA}

# recluster the central blob and get cluster markers
Rscript scripts/21_recluster_blob.R

# get gene universe annotation
grep -wf 13_recluster_blob/sycon_onlyBlob_geneUniverse.ls 09_gene_annotation/GOterms_OMA.tsv > 13_recluster_blob/sycon_onlyBlob_geneUniverse_GOterms.tsv

# get GO annotation for each cluster marker and the gene universe
for i in 13_recluster_blob/*/*ls; do grep -wf $i 09_gene_annotation/GOterms_OMA.tsv > "${i%%.*}"_GOterms.tsv; done

# perform GO enrichment for each cluster
for i in 13_recluster_blob/*/*_GOterms.tsv; do Rscript scripts/18_perform_GOenrich.R 13_recluster_blob/sycon_onlyBlob_geneUniverse_GOterms.tsv $i "${i%%.*}"_; done

# run hdWGCNA for reclusters
# TO COMPLETE

# run SAMap on re-clustered blob
python scripts/07_run_SAMap.py -s Scil,Aque -a pairwise -i 13_recluster_blob/04_SAMap/ -o 13_recluster_blob/04_SAMap/
python scripts/07_run_SAMap.py -s Scil,Slac -a pairwise -i 13_recluster_blob/04_SAMap/ -o 13_recluster_blob/04_SAMap/
python scripts/07_run_SAMap.py -s Scil,Aque,Slac -a pairwise -i 13_recluster_blob/04_SAMap/ -o 13_recluster_blob/04_SAMap/

# get SAMap statistics
for i in 13_recluster_blob/04_SAMap/*pkl; do python scripts/09_get_SAMap_mappingTables.py -p $i -o 13_recluster_blob/04_SAMap/01_mapping_scores -n 0; done
for i in 13_recluster_blob/04_SAMap/*pkl; do python scripts/10_get_SAMap_genePairs.py -p $i -o 13_recluster_blob/04_SAMap/02_gene_pairs -t 0.2; done
