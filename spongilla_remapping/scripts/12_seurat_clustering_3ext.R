#!/usr/bin/env Rscript

setwd("/data/evassvis/fn76/sycon/sycon_clusterAnnotation/spongilla_remapping")

library(Seurat)
library(SeuratExtend)
library(tidyverse)
library(ggvenn)


##########################
#     PRE-PROCESSING     #
##########################

# load 10X output
slac_3ext_remapped <- Read10X("02c_cellranger_count/03_aggregated/slac_3ext_aggr/outs/count/filtered_feature_bc_matrix/",
                              gene.column = 1) %>%
  CreateSeuratObject(project = "slac_3ext_remapped")
slac_3ext_remapped

# add a column in metadata for samples
slac_3ext_remapped[[]] <- slac_3ext_remapped[[]] %>%
  rownames_to_column(var = "cells") %>%
  separate(col = "cells", into = c("barcode", "sample"), sep = "-", remove = FALSE) %>%
  select(-barcode) %>%
  column_to_rownames(var = "cells") %>%
  mutate(sample = as_factor(sample))

# plot feature statistics before filtering
VlnPlot(slac_3ext_remapped,
        features = c("nFeature_RNA", "nCount_RNA"))
FeatureScatter(slac_3ext_remapped,
               feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "sample")

# filter based on features and counts
slac_3ext_remapped <- slac_3ext_remapped %>%
  subset(subset = nFeature_RNA > 200 & nFeature_RNA < 800 &
           nCount_RNA > 200 & nCount_RNA < 1500)

# plot feature statistics after filtering
feature_vln <- slac_3ext_remapped %>%
  VlnPlot(features = c("nFeature_RNA", "nCount_RNA"), group.by = "sample")
feature_scatter <- slac_3ext_remapped %>%
  FeatureScatter(feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "sample")

# create a panel
feature_panel <- ggpubr::ggarrange(feature_vln, feature_scatter)
feature_panel

# save the panel to file
ggsave("03_slac_remapped_clustering/02_plots/feature_plot_3ext_panel.pdf",
       feature_panel, device = cairo_pdf,
       width = 20, height = 8, dpi = 300, units = "in", bg = "white")
ggsave("03_slac_remapped_clustering/feature_plot_3ext_panel.png",
       feature_panel, device = "png",
       width = 20, height = 8, dpi = 300, units = "in", bg = "white")


######################
#     CLUSTERING     #
######################

slac_3ext_remapped <- slac_3ext_remapped %>%
  SCTransform(verbose = TRUE)

slac_3ext_remapped <- slac_3ext_remapped %>%
  RunPCA(verbose = TRUE)

ElbowPlot(slac_3ext_remapped, ndims = 50)

slac_3ext_remapped <- slac_3ext_remapped %>%
  RunUMAP(dims = 1:30, verbose = TRUE)

slac_3ext_remapped <- slac_3ext_remapped %>%
  FindNeighbors(dims = 1:30, verbose = TRUE) %>%
  FindClusters(verbose = TRUE, resolution = 2, cluster.name = "seurat_clusters_2") %>%
  FindClusters(verbose = TRUE, resolution = 1, cluster.name = "seurat_clusters_1")
  
Idents(object = slac_3ext_remapped) <- "seurat_clusters_2"
slac_3ext_remapped %>% saveRDS(file = "03_slac_remapped_clustering/slac_3ext_remapped_clustered.Rds")
slac_3ext_remapped <- readRDS(file = "03_slac_remapped_clustering/slac_3ext_remapped_clustered.Rds")

markers_slac_3ext <- FindAllMarkers(slac_3ext_remapped, group.by = "seurat_clusters", only.pos = TRUE)

markers_slac_3ext %>%
  filter(gene == "ENSLPGG00000010709")


###############################
#     SAVE RNA ASSAY ONLY     #
###############################

# this is to be used as input to SAMap
slac_3ext_remapped_RNA <- slac_3ext_remapped
DefaultAssay(slac_3ext_remapped_RNA) <- "RNA"
slac_3ext_remapped_RNA <- DietSeurat(slac_3ext_remapped_RNA, assays = "RNA")
scCustomize::as.anndata(x = slac_3ext_remapped_RNA, main_layer = "counts",
                        other_layers = NULL, file_path ="03_slac_remapped_clustering",
                        file_name = "slac_3ext_remapped_clustered.h5ad")


###################################
#     CHECK MARKER EXPRESSION     #
###################################

slac_3ext_remapped <- readRDS(file = "03_slac_remapped_clustering/slac_3ext_remapped_clustered.Rds")

Idents(object = slac_3ext_remapped) <- "seurat_clusters_2"

markers <- list("archaeocytes" = c("ENSLPGG00000008927", "ENSLPGG00000029352", "ENSLPGG00000010709", "ENSLPGG00000010195"),
                "choanocytes" = c("ENSLPGG00000014386"),
                "choanoblasts" = c("ENSLPGG00000009033"),
                "inc_pinacocytes" = c("ENSLPGG00000007262"))

SeuratExtend::DotPlot2(slac_3ext_remapped, #group.by = "seurat_clusters_2",
                       features = markers)

spongilla_3ext_allMarkers <- FindAllMarkers(slac_3ext_remapped,
                                            only.pos = TRUE)

top_markers_logFC1 <- spongilla_3ext_allMarkers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  # slice_head(n = 10) %>%
  ungroup()

top_markers_logFC1_list <- top_markers_logFC1 %>%
  group_by(cluster) %>%
  summarise(genes = list(unique(gene)), .groups = "drop")

common_markers <- crossing(cluster1 = top_markers_logFC1_list$cluster,
                           cluster2 = top_markers_logFC1_list$cluster) %>%
  filter(as.numeric(cluster1) < as.numeric(cluster2)) %>%
  rowwise() %>%
  mutate(
    genes1 = list(top_markers_logFC1_list$genes[top_markers_logFC1_list$cluster == cluster1][[1]]),
    genes2 = list(top_markers_logFC1_list$genes[top_markers_logFC1_list$cluster == cluster2][[1]]),
    shared_genes = list(intersect(genes1, genes2)),
    n_shared = length(shared_genes),
    n_unique_cluster1 = length(setdiff(genes1, genes2)),
    n_unique_cluster2 = length(setdiff(genes2, genes1)),
    n_shared_perc = n_shared / (n_unique_cluster1 + n_unique_cluster2 + n_shared) * 100,
  ) %>%
  select(cluster1, cluster2, n_shared, n_shared_perc, n_unique_cluster1, n_unique_cluster2, shared_genes)
common_markers

marker_tile <- common_markers %>%
  # filter(cluster1 %in% c("0", "1", "2", "3", "4", "5", "21", "23")) %>%
  # filter(cluster2 %in% c("0", "1", "2", "3", "4", "5", "21", "23")) %>%
  # # View()
  
  ggplot(aes(x = cluster2, y = cluster1, fill = n_shared)) +
  geom_tile(alpha = 0.5) +
  geom_text(aes(label = if_else(n_shared == 0, "", as.character(round(n_shared_perc,0)))),
            size = 2) +
  
  viridis::scale_fill_viridis() +
  scale_x_discrete(limits = rev) +
  
  labs(col = "Number of shared markers",
       title = "Shared markers between clusters\n(percentage within each tile)") +
  
  # scale_x_discrete(#limits = rownames(adj_mtx), 
  #                  labels = function(x) ifelse(x %in% photo_metacells, x, "")) +
  # scale_y_discrete(#limits = colnames(adj_mtx), 
  #                  labels = function(y) ifelse(y %in% photo_metacells, y, "")) +
  
  theme(aspect.ratio = 1,
        panel.background = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_text(size = 8))
marker_tile

ggsave(filename = "03_slac_remapped_clustering/02_plots/number_shared_markers_tile.png",
       marker_tile, device = "png",
       width = 8, height = 8, units = "in", dpi = 300)


ggvenn::ggvenn(list("cluster  15" = top_markers_logFC1 %>%
                      filter(cluster == 15) %>%
                      pull(gene),
                    "cluster  24" = top_markers_logFC1 %>%
                      filter(cluster == 24) %>%
                      pull(gene)),
               auto_scale = TRUE)


###################################
#     ASSIGN CLUSTER IDENTITY     #
###################################


cluster_identity <- c("0" = "Archaeocytes_1",
                      "1" = "Archaeocytes_2",
                      "2" = "Archaeocytes_3",
                      "3" = "Archaeocytes_4",
                      "4" = "Archaeocytes_5",
                      "5" = "Archaeocytes_4",
                      "6" = "Cluster_6",
                      "7" = "Cluster_7",
                      "8" = "Cluster_8",
                      "9" = "Myopeptidocytes",
                      "10" = "Cluster_10",
                      "11" = "Cluster_11",
                      "12" = "Cluster_12",
                      "13" = "Cluster_13",
                      "14" = "Metabolocytes",
                      "15" = "Myopeptidocytes",
                      "16" = "Cluster_16",
                      "17" = "Pinacocytes",
                      "18" = "Cluster_18",
                      "19" = "Cluster_19",
                      "20" = "Cluster_20",
                      "21" = "Archaeocytes_6",
                      "22" = "Cluster_22",
                      "23" = "Archaeocytes-like",
                      "24" = "Myopeptidocytes",
                      "25" = "Choanocytes/-blasts",
                      "26" = "Cluster_26",
                      "27" = "Cluster_27",
                      "28" = "Cluster_28",
                      "29" = "Cluster_29",
                      "30" = "Sclerocytes",
                      "31" = "Amoebocytes/Neuroid",
                      "32" = "Cluster_32",
                      "33" = "Cluster_33",
                      "34" = "Cluster_34",
                      "35" = "Basopinacocytes",
                      "36" = "Mesocytes",
                      "37" = "Granulocytes-like",
                      "38" = "Cluster_38",
                      "39" = "Cluster_39",
                      "40" = "Cluster_40",
                      "41" = "Mesocytes",
                      "42" = "Cluster_42")

slac_3ext_remapped[[]] <- slac_3ext_remapped[[]] %>%
  left_join(data.frame(cluster_identity) %>%
              rownames_to_column(var = "cluster_n"),
            by = join_by("seurat_clusters_2" == "cluster_n"))

umap <- DimPlot(slac_3ext_remapped, group.by = "cluster_identity") +
  theme_umap_arrows()
umap

umap@data <- umap@data %>%
  mutate(cluster_identity_reduced = str_replace(cluster_identity, "_[0-9]+$", ""))

umap_cell_types <- umap@data %>%
  ggplot(aes(umap_1, umap_2)) +
  
  geom_point(data = subset(umap@data, grepl("Cluster", cluster_identity_reduced)),
             colour = "grey85", size = 1) +
  geom_point(data = subset(umap@data, !grepl("Cluster", cluster_identity_reduced)),
             aes(colour = as.factor(cluster_identity_reduced)),
             size = 1) +
  
  # scale_colour_discrete(name = "",
  labs(colour = "Main cell types") +
  
  theme_void() +
  theme_umap_arrows()
umap_cell_types


umap_cell_ID <- DimPlot2(slac_3ext_remapped, features = "seurat_clusters_2", pt.size = 1,
         label = TRUE, repel = TRUE, box = TRUE, label.color = "black",
         theme = list(labs(title = ""), NoLegend(), theme_umap_arrows()))

umap_cell_ID

panel <- ggpubr::ggarrange(umap_cell_ID, umap_cell_types, align = "hv")

ggsave(filename = "03_slac_remapped_clustering/02_plots/panel_umaps_cell_atlas.png",
       panel, device = "png",
       width = 20, height = 8, units = "in", dpi = 300)

