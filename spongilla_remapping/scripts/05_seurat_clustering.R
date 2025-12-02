#!/usr/bin/env Rscript

setwd("/data/evassvis/fn76/sycon/sycon_clusterAnnotation/spongilla_remapping")

library(Seurat)
library(SeuratExtend)
library(tidyverse)


##########################
#     PRE-PROCESSING     #
##########################

# load 10X output
slac_remapped <- Read10X("02a_cellranger_count/03_aggregated/slac_aggr/outs/count/filtered_feature_bc_matrix/",
                         gene.column = 1) %>%
  CreateSeuratObject(project = "slac_remapped")
slac_remapped

# add a column in metadata for samples
slac_remapped[[]] <- slac_remapped[[]] %>%
  rownames_to_column(var = "cells") %>%
  separate(col = "cells", into = c("barcode", "sample"), sep = "-", remove = FALSE) %>%
  select(-barcode) %>%
  column_to_rownames(var = "cells") %>%
  mutate(sample = as_factor(sample))

# plot feature statistics before filtering
VlnPlot(slac_remapped,
        features = c("nFeature_RNA", "nCount_RNA"), group.by = "sample")
FeatureScatter(slac_remapped,
               feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "sample")

# filter based on features and counts
slac_remapped <- slac_remapped %>%
  subset(subset = nFeature_RNA > 200 & nFeature_RNA < 800 &
           nCount_RNA > 200 & nCount_RNA < 1500)

# plot feature statistics after filtering
feature_vln <- slac_remapped %>%
  VlnPlot(features = c("nFeature_RNA", "nCount_RNA"), group.by = "sample")
feature_scatter <- slac_remapped %>%
  FeatureScatter(feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "sample")

# create a panel
feature_panel <- ggpubr::ggarrange(feature_vln, feature_scatter)
feature_panel

# save the panel to file
ggsave("03_slac_remapped_clustering/feature_plot_panel.pdf",
       feature_panel, device = cairo_pdf,
       width = 20, height = 8, dpi = 300, units = "in", bg = "white")
ggsave("03_slac_remapped_clustering/feature_plot_panel.png",
       feature_panel, device = "png",
       width = 20, height = 8, dpi = 300, units = "in", bg = "white")


######################
#     CLUSTERING     #
######################

slac_remapped <- slac_remapped %>%
  SCTransform(verbose = TRUE)

slac_remapped <- slac_remapped %>%
  RunPCA(verbose = TRUE)

ElbowPlot(slac_remapped, ndims = 50)

slac_remapped <- slac_remapped %>%
  RunUMAP(dims = 1:30, verbose = TRUE)

slac_remapped <- slac_remapped %>%
  FindNeighbors(dims = 1:30, verbose = TRUE) %>%
  FindClusters(verbose = TRUE, resolution = 1)

slac_remapped %>% saveRDS(file = "03_slac_remapped_clustering/slac_remapped_clustered.Rds")
slac_remapped <- readRDS(file = "03_slac_remapped_clustering/slac_remapped_clustered.Rds")

DimPlot(slac_remapped)


###############################
#     SAVE RNA ASSAY ONLY     #
###############################

# this is to be used as input to SAMap
slac_remapped_RNA <- slac_remapped
DefaultAssay(slac_remapped_RNA) <- "RNA"
slac_remapped_RNA <- DietSeurat(slac_remapped_RNA, assays = "RNA")
scCustomize::as.anndata(x = slac_remapped_RNA, main_layer = "counts",
                        other_layers = NULL, file_path ="03_slac_remapped_clustering",
                        file_name = "slac_remapped_clustered.h5ad")


############################################
#     TRY TO ANNOTATE BASED ON MARKERS     #
############################################

slac_original <- slac_original %>%
  SCTransform(verbose = TRUE)

slac_original <- slac_original %>%
  RunPCA(verbose = TRUE)

ElbowPlot(slac_original, ndims = 50)

slac_original <- slac_original %>%
  RunUMAP(dims = 1:30, verbose = TRUE)

# slac_original <- slac_original %>%
#   FindNeighbors(dims = 1:30, verbose = TRUE) %>%
#   FindClusters(verbose = TRUE, resolution = 1)

saveRDS(slac_original, "03_slac_remapped_clustering/slac_original_clustered.Rds")
slac_original <- readRDS("03_slac_remapped_clustering/slac_original_clustered.Rds")

slac_original[[]] <- slac_original[[]] %>%
  mutate(clusterID = as_factor(clusterID))
SeuratExtend::DimPlot2(slac_original, group.by = "cell_type_abbreviation", repel = TRUE,
                       label = TRUE, box = TRUE, theme = theme_umap_arrows())

markers_original <- FindAllMarkers(slac_original, group.by = "cell_type_abbreviation", only.pos = TRUE)

markers_remapped <- FindAllMarkers(slac_remapped, group.by = "seurat_clusters", only.pos = TRUE)
colplot <- markers_remapped %>%
  group_by(cluster) %>%
  arrange(-avg_log2FC) %>%
  slice_head(n = 50) %>%
  ungroup() %>%
  left_join(read.table("03_slac_remapped_clustering/01_diamond_transcriptome_genome/slac_genomeVStranscriptome.tsv",
                       header = FALSE, sep = "\t", quote = "", na.strings = "") %>%
              group_by(V1) %>%
              slice_max(V3) %>%
              ungroup() %>%
              select(V1, V2) %>%
              mutate(V1 = str_replace_all(V1, "ENSLPGP", "ENSLPGG"),
                     V1 = str_replace_all(V1, "\\.1$", "")),
            by = join_by("gene" == "V1"),
            relationship = "many-to-many") %>%
  left_join(markers_original %>%
              select(gene, cluster, p_val_adj, avg_log2FC) %>%
              filter(p_val_adj <= .05),
            by = join_by("V2" == "gene"),
            relationship = "many-to-many") %>%
  drop_na() %>%
  
  group_by(cluster.x, cluster.y) %>%
  count() %>%
  
  ggplot(aes(x = cluster.x, y = n, fill = cluster.y)) +
    geom_col() +
  theme_minimal()

colplot

ggsave("04_SAMap/slac_integrated_clusterComposition.pdf",
       colplot, device = cairo_pdf,
       width = 10, height = 10, dpi = 300, units = "in", bg = "white")

  
markers_original %>%
  filter(p_val_adj < .05) %>%
  filter(gene %in% c("c95872-g2","c39492-g1","c104944-g2","c104138-g2","c95050-g1","c103117-g1","c104445-g4"))
