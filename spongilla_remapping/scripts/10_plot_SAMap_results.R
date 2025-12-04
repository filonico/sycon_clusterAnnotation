#!/usr/bin/env Rscript

setwd("/data/evassvis/fn76/sycon/sycon_clusterAnnotation/spongilla_remapping")

library(Seurat)
library(SeuratExtend)
library(tidyverse)
library(ggsankeyfier)
library(ggpubr)
library(RColorBrewer)
library(purrr)
library(igraph)


#####################
#     UMAP SAMAP    #
#####################

samap <- sceasy::convertFormat("04_SAMap/slac_integrated_originalVSremapped.h5ad",
                               from = "anndata", to = "seurat")

umap <- DimPlot(samap, group.by = "slacOriginal_cell_type", pt.size = 0.2) +
  # scale_color_manual(labels = c("mapped to\ntranscriptome", "mapped to\ngenome"),
  #                    values = c("#ff0f7b", "#f89b29")) +
  ggtitle("Integrated\nSpongilla datasets") +
  theme_umap_arrows()
umap

umap <- umap@data %>%
  ggplot(aes(UMAP_1, UMAP_2)) +
  
  geom_point(data = subset(umap@data, grepl("unassigned", slacOriginal_cell_type)),
             colour = "grey80", size = 1) +
  geom_point(data = subset(umap@data, !grepl("unassigned", slacOriginal_cell_type)),
             aes(colour = as.factor(slacOriginal_cell_type)),
             size = 1) +

  # scale_colour_discrete(name = "",
  labs(title = "Integrated Spongilla dataset", colour = "Cell types in the\noriginal dataset") +
  
  theme_void() +
  theme_umap_arrows()
umap

SeuratExtend::DimPlot2(samap, features = c("slacOriginal-c92270-g1", "slacRemapped-ENSLPGG00000029352"))

ggsave(filename = "04_SAMap/slac_integrated_originalVSremapped_UMAP_remappedGrey.pdf",
       umap, device = cairo_pdf,
       width = 12, height = 8, units = "in", dpi = 300, bg = "white")

ggsave(filename = "04_SAMap/slac_integrated_originalVSremapped_UMAP_remappedGrey.png",
       umap, device = "png",
       width = 12, height = 8, units = "in", dpi = 300, bg = "white")


###################
#     HEATMAP     #
###################

heatmap_mappingScores <- read.table("04_SAMap/03_SAMap_statistics/originalVSremapped_100topCells_samapMappingTable.tsv", sep = "\t",
           header = TRUE) %>%
  pivot_longer(-X, names_to = "to", values_to = "weight") %>%
  filter(grepl("slacOriginal", X) &
           grepl("slacRemapped", to)) %>%
  mutate(weight = if_else(weight < 0.2, 0, weight)) %>%
  
  ggplot(aes(x = X, y = to, fill = weight)) +
  geom_tile(col = "white", alpha = 0.5) +
  geom_text(aes(label = if_else(weight == 0, "", as.character(round(weight, 1)))),
            size = 2) +
  
  
  viridis::scale_fill_viridis() +
  scale_x_discrete(limits = rev) +

  labs(fill = "Mappin score",
       title = "SAMap mapping scores",
       x = "Original clusters", y = "Clusters on the remapped dataset") +
  
  theme(aspect.ratio = 1,
        panel.background = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8))

heatmap_mappingScores

ggsave(filename = "04_SAMap/slac_integrated_originalVSremapped_mappingScores_tile.png",
       heatmap_mappingScores, device = "png",
       width = 9, height = 9, units = "in", dpi = 300, bg = "white")
ggsave(filename = "04_SAMap/slac_integrated_originalVSremapped_mappingScores_tile.pdf",
       heatmap_mappingScores, device = cairo_pdf,
       width = 9, height = 9, units = "in", dpi = 300, bg = "white")

