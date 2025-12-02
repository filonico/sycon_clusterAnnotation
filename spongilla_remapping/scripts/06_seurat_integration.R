#!/usr/bin/env Rscript

setwd("/data/evassvis/fn76/sycon/sycon_clusterAnnotation/spongilla_remapping")

library(Seurat)
library(SeuratExtend)
library(tidyverse)


####################################
#     INEGRATION WITH ORIGINAL     #
####################################

# load original slac dataset
slac_original <- sceasy::convertFormat("../sycon_cluster_analyses/04_preprocessed_scRNAseqs/Slac_cellFiltered.h5ad",
                                       from = "anndata", to = "seurat")

# load diamond best hits between genome and transcriptome proteins
diamond_best_hits_tVSg <- read.table("03_slac_remapped_clustering/01_diamond_transcriptome_genome/slac_transcriptomeVSgenome.tsv",
                                     header = FALSE, sep = "\t", quote = "", na.strings = "") %>%
  group_by(V1) %>%
  slice_max(V3) %>%
  slice_min(V11) %>%
  slice_max(V4) %>%
  ungroup() %>%
  filter(V3 >= 85) 
  # select(c(V1, V2)) %>%
  # rename("transcriptome_feature_name" = "V1", "genome_feature_name" = "V2")

diamond_best_hits_gVSt <- read.table("03_slac_remapped_clustering/01_diamond_transcriptome_genome/slac_genomeVStranscriptome.tsv",
                                     header = FALSE, sep = "\t", quote = "", na.strings = "") %>%
  group_by(V1) %>%
  slice_max(V3) %>%
  slice_min(V11) %>%
  slice_max(V4) %>%
  ungroup() %>%
  filter(V3 >= 85)

diamond_reciprocal_bestHits <- diamond_best_hits_tVSg %>%
  left_join(diamond_best_hits_gVSt, by = c("V2" = "V1"), suffix = c("tVSg", "gVSt"),
            relationship = "many-to-many") %>%
  filter(V1 == V2gVSt) %>%
  # filter(V2 == "ENSLPGP00000022701.1") %>%
  # View()
  # select(V2) %>%
  group_by(V2) %>% 
  filter(n()==1) %>%
  select(V1, V2) %>%
  ungroup()

# recode original slac dataset to have the genome gene names
slac_original_mtx <- as.matrix(slac_original@assays$RNA@counts) %>%
  as.data.frame() %>%
  rownames_to_column(var = "features") %>%
  left_join(diamond_reciprocal_bestHits, by = c("features" = "V1")) %>%
  mutate(genome_feature_name = coalesce(V2, features),
         genome_feature_name = str_replace_all(genome_feature_name, "\\.1$", ""),
         genome_feature_name = str_replace_all(genome_feature_name, "ENSLPGP", "ENSLPGG")) %>%
  select(-c(features,V2)) %>%
  column_to_rownames(var = "genome_feature_name")

slac_original_recoded <- CreateSeuratObject(counts = slac_original_mtx,
                                            project = "slac_original_recoded",
                                            min.cells = 0, min.features = 0) %>%
  AddMetaData(slac_original[[c("clusterID", "cell_type", "cell_type_abbreviation", "cell_type_family", "cell_type_newName")]])
slac_original_recoded

objects_to_integrate <- list(slac_original_recoded,
                             slac_remapped_RNA)

# normalize and identify variable features for each dataset independently
objects_to_integrate <- lapply(X = objects_to_integrate, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 5000)
})

# select features that are repeatedly variable across datasets for integration
repeatedly_variableFeatures <- SelectIntegrationFeatures(object.list = objects_to_integrate,
                                                         nfeatures = 2000)

# run PCA on each dataset using repeatedly variable features
objects_to_integrate <- lapply(X = objects_to_integrate, FUN = function(x) {
  x <- ScaleData(x, features = repeatedly_variableFeatures, verbose = FALSE)
  x <- RunPCA(x, features = repeatedly_variableFeatures, verbose = FALSE)
})

# integrate datasets
integration_anchors <- FindIntegrationAnchors(object.list = objects_to_integrate,
                                              anchor.features = repeatedly_variableFeatures,
                                              reduction = "rpca", k.anchor = 25)

saveRDS(integration_anchors, file = "03_slac_remapped_clustering/integration_anchors.Rds")

# this command creates an 'integrated' data assay
integrated_dataset <- IntegrateData(anchorset = integration_anchors)

saveRDS(integrated_dataset, file = "03_slac_remapped_clustering/integrated_dataset.Rds")

# Run the standard workflow for visualization and clustering
integrated_dataset <- ScaleData(integrated_dataset, verbose = FALSE)
integrated_dataset <- RunPCA(integrated_dataset)
ElbowPlot(integrated_dataset, ndims = 50)
integrated_dataset <- RunUMAP(integrated_dataset, reduction = "pca", dims = 1:30)
integrated_dataset <- FindNeighbors(integrated_dataset, reduction = "pca", dims = 1:30)
integrated_dataset <- FindClusters(integrated_dataset)
saveRDS(integrated_dataset, file = "03_slac_remapped_clustering/integrated_dataset_UMAP.Rds")
SeuratExtend::DimPlot2(integrated_dataset, split.by = "orig.ident")
