#!/bin/bash


##########################################
#     TIDY UP N. VECTENSIS PROTEOMES     #
##########################################

# two proteomes were merged, as in the UMI table there were two kinds of features IDs: one proteome was taken from Matt, the other from ENSEMBL
# merge the two Nvec proteomes and let headers be the same as in UMI tables
zcat 00_input/genomic_data/N_vectensis/Nematostella_vectensis.ASM20922v1.pep.all.fa.gz 00_input/genomic_data/N_vectensis/Nvec_Technau.protein.fasta.gz | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | tail -n +2 | sed -E 's/CnidariaNematostellavectensis_//g; s/^.+NEMVEDRAFT_/>/g; s/ .+$//; s/^>/>Nvec-/' > 00_input/genomic_data/N_vectensis/merged_proteome.faa

# the file Nvec_features.ls was created extracted features from the Rds file using R. Find the command below
# readRDS("00_input/scRNAseqs_rawTables/N_vectensis/Nvec_adult_sc_UMI_counts.RDS") %>% CreateSeuratObject(min.cells = 0, min.features = 0) %>% Features() %>% write.table(file = "/alice-home/2/f/fn76/ANALYSIS/sycon_clusterAnnotation/00_input/scRNAseqs_rawTables/N_vectensis/Nvec_features.ls", sep = "", row.names = FALSE, quote = FALSE)

# downsample the merged proteome with genes present in the UMI tables
python scripts/02_extract_sequences_from_fasta.py -l 00_input/scRNAseqs_rawTables/N_vectensis/Nvec_features.ls -f 00_input/genomic_data/N_vectensis/merged_proteome.faa -o 01_proteomes/Nvec.faa

mv 01_proteomes/Nvec_notFound.ls 00_input/genomic_data/N_vectensis/


#####################################
#     TIDY UP S. LACUSTRIS DATA     #
#####################################

# proteome taken from Matt
# remove annotation from gene IDs in UMI matrix
zcat 00_input/scRNAseqs_rawTables/S_lacustris/GSE134912_Slac_spongilla_10x_count_matrix.txt.gz | sed -E 's/ [^\t]+//' | gzip -v -9 - > 00_input/scRNAseqs_rawTables/S_lacustris/GSE134912_Slac_spongilla_10x_count_matrix_edited.txt.gz

# edit gene IDs in the proteome fasta
zcat 00_input/genomic_data/S_lacustris/Spongilla.fasta.gz | sed -E 's/-/_/' > 01_proteomes/Slac.faa


#############################################
#     TIDY UP A. QUEENSLANDICA PROTEOME     #
#############################################

# proteome taken from Matt
# edit gene IDs in the proteome fasta
zcat 00_input/genomic_data/A_queenslandica/Aque.fasta.gz | sed -E 's/-/_/' > 01_proteomes/Aque.faa


#############################################
#     TIDY UP C. HEMISPHAERICA PROTEOME     #
#############################################

#### REMOVED ####

# proteome downloaded from https://data.caltech.edu/records/6vy5w-0a884 (21/02/2025)
# edit gene IDs in the proteome fasta
zcat 00_input/genomic_data/C_hemisphaerica/20201030_cdhit_95.fasta.transdecoder.pep.gz | sed -E 's/ .+$//' > 01_proteomes/Chem.faa


#########################################
#     TIDY UP T. ADHAERENS PROTEOME     #
#########################################

# proteome downloaded form https://github.com/xgrau/placozoa-cell-type-evolution-code/tree/master (21/02/2025)
# edit gene IDs in the proteome fasta
sed -E 's/TriadT/TriadG/; s/\..+$//' placozoa-cell-type-evolution-code/data/reference/Tadh_long.pep.fasta > 01_proteomes/Tadh.faa


##########################################
#     TIDY UP S. PISTILLATA PROTEOME     #
##########################################

# the proteomes were downloaded from NBCI (GCF_002571385.2) and from ReefGenomics (24/02/2025)
# merge the two proteomes and make the headers match those of the UMI table
zcat 00_input/genomic_data/S_pistillata/*gz | sed -E 's/ .+$//; s/\.1/_1/' | sed -E '/^>Spis/ s/$/_1/' | sed 's/>XP/>Spis_XP/; s/_/-/g' > 00_input/genomic_data/S_pistillata/Spis_merged_proteins_renamed.faa

# the file Spis_features.ls was created extracted features from the Rds file using R. Find the command below
# readRDS("00_input/scRNAseqs_rawTables/S_pistillata/Spis_adult_sc_UMI_counts.RDS") %>% CreateSeuratObject(min.cells = 0, min.features = 0) %>% Features() %>% write.table(file = "/alice-home/2/f/fn76/ANALYSIS/sycon_clusterAnnotation/00_input/scRNAseqs_rawTables/S_pistillata/Spis_features.ls", sep = "", row.names = FALSE, quote = FALSE)

# extract genes matching the list of features
python scripts/02_extract_sequences_from_fasta.py -l 00_input/scRNAseqs_rawTables/S_pistillata/Spis_features.ls -f 00_input/genomic_data/S_pistillata/Spis_merged_proteins_renamed.faa -o 01_proteomes/Spis.faa

mv 01_proteomes/Spis_notFound.ls 00_input/genomic_data/S_pistillata/


########################################
#     TIDY UP H. VULGARIS PROTEOME     #
########################################

# the proteome was downloaded from https://research.nhgri.nih.gov/hydra/ (24/02/2025)
# merge the two proteomes and make the headers match those of the UMI table
zcat 00_input/genomic_data/H_vulgaris/hydra2.0_genemodels.aa.gz | grep --no-group-separator -A1 t1 | sed -E 's/^>.+\.g/>Hvul-g/; s/\..+$/-1/' > 00_input/genomic_data/H_vulgaris/hydra2.0_genemodels_renamed.aa

# the file Hvul_features.ls was created extracted features from the Rds file using R. Find the command below
# readRDS("00_input/scRNAseqs_rawTables/H_vulgaris/Hvul_sc_UMI_counts.RDS") %>% CreateSeuratObject(min.cells = 0, min.features = 0) %>% Features() %>% write.table(file = "/alice-home/2/f/fn76/ANALYSIS/sycon_clusterAnnotation/00_input/scRNAseqs_rawTables/H_vulgaris/Hvul_features.ls", sep = "", row.names = FALSE, quote = FALSE)

# extract genes matching the list of features
python scripts/02_extract_sequences_from_fasta.py -l 00_input/scRNAseqs_rawTables/H_vulgaris/Hvul_features.ls -f 00_input/genomic_data/H_vulgaris/hydra2.0_genemodels_renamed.aa -o 01_proteomes/Hvul.faa

mv 01_proteomes/Hvul_notFound.ls 00_input/genomic_data/H_vulgaris/


##########################
#     TIDY UP X. SP.     #
##########################

# the proteome was downloaded from https://cmo.carnegiescience.edu/endosymbiosis/genome/ (24/02/2025)
zcat 00_input/genomic_data/X_sp/xenSp1.proteins.fa.gz | grep --no-group-separator -A1 T1 | sed -E 's/-T1.+$//; s/>Xe_/>Xesp-/' > 00_input/genomic_data/X_sp/xenSp1.proteins.renamed.fa

# the file Xesp_features.ls was created extracted features from the Rds file using R. Find the command below
# readRDS("00_input/scRNAseqs_rawTables/X_sp/Xesp_sc_UMI_counts.RDS") %>% CreateSeuratObject(min.cells = 0, min.features = 0) %>% Features() %>% write.table(file = "/alice-home/2/f/fn76/ANALYSIS/sycon_clusterAnnotation/00_input/scRNAseqs_rawTables/X_sp/Xesp_features.ls", sep = "", row.names = FALSE, quote = FALSE)

# extract genes matching the list of features
python scripts/02_extract_sequences_from_fasta.py -l 00_input/scRNAseqs_rawTables/X_sp/Xesp_features.ls -f 00_input/genomic_data/X_sp/xenSp1.proteins.renamed.fa -o 01_proteomes/Xesp.faa

mv 01_proteomes/Xesp_notFound.ls 00_input/genomic_data/X_sp/



##########################################
#     ONELINE AND REMOVE UNDERSCORES     #
##########################################

# underscores are replaced with "-" because Seurat does not like them
for i in 01_proteomes/*faa; do FILENAME="${i%.*}"; FILEEXTENSION="${i##*.}"; awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $i | tail -n +2 | sed -E 's/_/-/g' > "$FILENAME"_ol."$FILEEXTENSION"; done
