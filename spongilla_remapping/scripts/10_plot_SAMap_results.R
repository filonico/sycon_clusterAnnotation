#!/usr/bin/env Rscript

setwd("/data/evassvis/fn76/sycon/sycon_clusterAnnotation/spongilla_remapping")

library(Seurat)
library(SeuratExtend)
library(tidyverse)
library(ggsankeyfier)
library(ggpubr)
library(RColorBrewer)
library(purrr)


#####################
#     FUNCTIONS     #
#####################

# function to prepare the MappingScore SAMap matrix for plottin
tidyup_dataframe <- function(filename, threshold) {
  
  # define the discrete categories of mapping scores
  value_intervals <- c(0, 0.25, 0.5, 0.75, 1)
  value_interval_labels <- c("\u2264 0.25", "\u2264 0.50", "\u2264 0.75", "\u2264 1.00")
  
  dataframe <- read.table(filename, header = TRUE, sep = "\t") %>%
    
    # create a long-format dataframe
    pivot_longer(-X, names_to = "target", values_to = "value") %>%
    rename(source = X) %>%
    
    # filter out rows with mapping scores below the threshold
    filter(value >= threshold) %>%
    
    # split the column of cell cluster IDs into two:
    # one with the species name and another with the cell cluster name
    separate(source, into = c("x", "group"), sep = "_", extra = "merge") %>%
    separate(target, into = c("next_x", "next_group"), sep = "_", extra = "merge") %>%
    
    # create a column with unique edge IDs (necessary for the sankey plot function)
    mutate(edge_id = paste(pmin(group, next_group), pmax(group, next_group), sep = "_")) %>%
    
    # tidy up cell cluster names
    mutate(edge_id = str_replace_all(edge_id, "_", "."),
           edge_id = str_replace_all(edge_id, "-", "."),
           edge_id = str_replace_all(edge_id, " ", "."),
           group = str_replace_all(group, "_", " "),
           next_group = str_replace_all(next_group, "_", " ")) %>%
    
    # create a column with discrete categories for mapping scores
    mutate(category = as_factor(cut(value, breaks = value_intervals,
                                    labels = as.factor(value_interval_labels),
                                    right = FALSE))) %>%
    
    # create a column with cleaned cell type names, in order to color code nodes
    mutate(cell_type_col = str_remove(group, "\\s.*"),
           category = factor(category, levels = c("\u2264 0.25", "\u2264 0.50", "\u2264 0.75", "\u2264 1.00")))
  
  return(dataframe)
}

# function to produce a sankey plot out of the tidied-up dataframe
plot_sankey <- function(dataframe, species, species_colors) {
  
  # define element and text position in the plot
  pos <- position_sankey(align = "center", v_space = "auto")
  pos_text_left <- position_sankey(align = "center", v_space = "auto", nudge_x = -0.09)
  pos_text_right <- position_sankey(align = "center", v_space = "auto", nudge_x = 0.09)
  
  colourCount <- length(unique(dataframe$group))
  getPalette <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
  
  dataframe <- dataframe %>%
    mutate(connector = case_when(
      x == "Transcriptome" ~ "from",              # If x is Sycon, connector is "from"
      # x < next_x ~ "from",                # Default condition: x is less than next_x
      TRUE ~ "to"                         # Otherwise, connector is "to"
    ),
    
    # sort species names
    x_order = case_when(
      x == "Transcriptome" ~ 1,          # If x is "Sycon", give it the first priority
      TRUE ~ 3                    # Otherwise, assign the third priority (for all other species)
    ),
    
    x = factor(x, levels = species)) %>%
    
    arrange(x_order, x)
  
  number_species1 <- sum(dataframe$x == species[1])
  
  plot <- dataframe %>%
    
    ggplot(aes(x = x, y = value, group = group,
               connector = connector, edge_id = edge_id)) +
    
    # add sankey edges, colored by categories
    ggsankeyfier::geom_sankeyedge(aes(fill = category),
                                  position = pos, alpha = 0.7) +
    scale_fill_grey(name = "Mapping\nscores",
                    start = 0.9, end = 0,
                    guide = guide_legend(reverse = TRUE),
                    drop = FALSE
    ) +
    
    # add a new scale fill layer
    ggnewscale::new_scale_fill() +
    
    # add edge annotation with cell cluster names
    # left labels
    geom_text(aes(label = c(group[1:number_species1], rep(NA, length(group)-number_species1))), hjust = 1,
              fontface = "bold",
              col = "grey15", stat = "sankeynode", position = pos_text_left, cex = 3) +
    # right labels
    geom_text(aes(label = c(rep(NA, length(group)-number_species1), group[(number_species1+1):length(group)])), hjust = 0,
              fontface = "bold",
              col = "grey15", stat = "sankeynode", position = pos_text_right, cex = 3) +
    
    # add sankey nodes colored by cell type
    ggsankeyfier::geom_sankeynode(aes(#col = as.factor(group),
      fill = as.factor(cell_type_col),
      # fill = stage(as.factor(group),
      #              after_scale = alpha(fill, 0.5)))
    ),
    col = "grey15",
    linewidth = 0.4, position = pos) +
    scale_fill_manual(values = getPalette(colourCount),
                      #values = species_colors,
                      guide = "none", aesthetics = c("color", "fill")) +
    
    # prevent clipping of plot elements
    coord_cartesian(clip = 'off') +
    
    # define plot themes
    theme_minimal() +
    theme(axis.text.x = element_text(color = species_colors,
                                     face = "italic",
                                     size = 12),
          axis.text.y = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank(),
          legend.title = element_text(hjust = 0.5, face = "bold"))
  
  return(plot) 
}

fromTable_toSankey <- function(file_path, prefix_sp1, prefix_sp2, name_sp1, name_sp2, threshold = 0.2) {
  
  df <- tidyup_dataframe(file_path, threshold) %>%
    mutate(
      x = str_replace_all(x, prefix_sp1, name_sp1),
      x = str_replace_all(x, prefix_sp2, name_sp2),
      next_x = str_replace_all(next_x, prefix_sp1, name_sp1),
      next_x = str_replace_all(next_x, prefix_sp2, name_sp2)
    )
  
  plot_sankey(df, c(name_sp1, name_sp2), rep("grey15", 2))
  
}


#####################
#     UMAP SAMAP    #
#####################

HH23HhonTadhTrH2_samap <- sceasy::convertFormat("06_SAMap/HH23HhonTadhTrH2_leiden3Clusters_samap.h5ad",
                                                from = "anndata", to = "seurat")
continous_scale <- c(low = "lightblue", high = "darkblue")

HH23HhonTadhTrH2_samap[[]] <- HH23HhonTadhTrH2_samap[[]] %>%
  rownames_to_column(var = "cells") %>%
  pivot_longer(cols = ends_with("_cell_type"),
               names_to = "species_col",
               values_to = "value") %>%
  # keep only the matching species-column pairs
  filter(species_col == paste0(species, "_cell_type")) %>%
  mutate(cell_type = ifelse(value == "unassigned" | is.na(value),
                               paste0(species, "_unassigned"),
                               value)) %>%
  select(-species_col, -value, -ends_with("_Metacell_ale")) %>%
  
  column_to_rownames(var = "cells")


SeuratExtend::DimPlot2(HH23HhonTadhTrH2_samap,
                       features = c("Hhon-Hhon-g05380.t1","Hhon-Hhon-g05381.t1",
                                    "HH23-HoiH23-PlH23-000307","HH23-HoiH23-PlH23-000306",
                                    "Tadh-Tadh-TriadG26040","Tadh-Tadh-TriadG26009",
                                    "TrH2-TrH2-TrispH2-008979","TrH2-TrH2-TrispH2-008980"),
                       pt.size = 0.1, nrow = 4, cols = c(continuous = continous_scale)) +
  theme_umap_arrows()
umap

umap <- umap@data %>%
  ggplot(aes(UMAP_1, UMAP_2)) +
  
  geom_point(data = subset(umap@data, grepl("unassigned", slacOriginal_cell_type)),
             colour = "grey80",
             size = 1, alpha = 0.6) +
  geom_point(data = subset(umap@data, !grepl("unassigned", slacOriginal_cell_type)),
             aes(colour = as.factor(slacOriginal_cell_type)),
             size = 1, alpha = 0.6) +
  
  # scale_colour_discrete(name = "",
                        # labels = c("Mapped to genome", "Mapped to transcriptome")) +
  
  ggtitle("Integrated Spongilla datasets") +
  
  theme_bw()

umap

ggsave(filename = "04_SAMap/slac_integrated_originalVSremapped_UMAP_remappedGrey.pdf",
       umap, device = cairo_pdf,
       width = 12, height = 8, units = "in", dpi = 300)

ggsave(filename = "04_SAMap/slac_integrated_originalVSremapped_UMAP_remappedGrey.png",
       umap, device = "png",
       width = 12, height = 8, units = "in", dpi = 300)


########################
#     SANKEY PLOTS     #
########################

fromTable_toSankey(
  "04_SAMap/03_SAMap_statistics/originalVSremapped_100topCells_samapMappingTable.tsv",
  "slacOriginal", "slacRemapped",
  "Transcriptome", "Genome", 0.3
)
