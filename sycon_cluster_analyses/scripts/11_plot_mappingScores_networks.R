setwd("/lustre/alice3/data/evassvis/fn76/sycon/sycon_clusterAnnotation/sycon_cluster_analyses")

library(igraph)
library(ggraph)
library(graphlayouts)
library(ggforce)
library(tidyverse)
library(ggrepel)
library(pals)


#####################
#     FUNCTIONS     #
#####################

# function to plot the network with metacells
get_graph_mtx_metacell <- function(file, threshold) {
  adj_mtx <- read.table(file, header = TRUE) %>%
    as.matrix()
  
  # set to 0 all the values below a threshold
  adj_mtx[adj_mtx < threshold] <- 0
  
  # define all-0 rows and columns (i.e., non-mapping metacells)
  zero_rows <- rowSums(adj_mtx) == 0
  zero_cols <- colSums(adj_mtx) == 0
  
  # remove rows/columns
  adj_mtx_reduced <- adj_mtx[!zero_rows, !zero_cols, drop = FALSE] # remove all 0-columns/rows, regardless
  
  # build a graph from matrix
  graph_mtx <- adj_mtx_reduced %>%
    igraph::graph_from_adjacency_matrix(weighted = TRUE, mode = "undirected", diag = FALSE)
  
  return(graph_mtx)  
}

plot_network_metacell <- function(graph_mtx) {
  
  # define colour per species
  V(graph_mtx)$species <- substr(V(graph_mtx)$name, 1, 4)
  species_fills <- c(
    "Slac" = "#b5c9da",
    "Aque" = "#c0a3c6",
    "Scil" = "#f2bba4"
  )
  
  species_cols <- c(
    "Slac" = "#07589d",
    "Aque" = "#722376",
    "Scil" = "#dd402c"
  )
  
  # Bin edge weights into categories
  value_intervals <- c(0, 0.5, 0.75, 1)
  value_interval_labels <- c("\u2264 0.50", "\u2264 0.75", "\u2264 1.00")
  
  E(graph_mtx)$weight_bin <- cut(
    E(graph_mtx)$weight,
    breaks = value_intervals,
    labels = value_interval_labels,
    include.lowest = TRUE,
    right = TRUE
  )
  
  # plot the graph
  graph_plot <- graph_mtx %>%
    ggraph::ggraph(layout = "stress") +
    
    # define edge style
    geom_edge_link0(aes(edge_width = weight_bin, edge_colour = weight_bin),
                    edge_alpha = 0.7) +
    
    # define node styles
    geom_node_point(aes(fill = species, col = species), size = 5, shape = 21) +
    
    # define text style (display only that of photo-metacells)
    geom_node_text(aes(label = name),
                   size = 3, repel = TRUE, max.overlaps = Inf) +
    
    scale_edge_colour_grey(name = "Mapping\nscores", start = 0.9, end = 0,
                           guide = guide_legend(reverse = TRUE), drop = FALSE) +
    
    scale_edge_width_manual(name = "Mapping\nscores", guide = guide_legend(reverse = TRUE),
                            values = c(
                              # "\u2264 0.25" = 0.5,
                              "\u2264 0.50" = 1,
                              "\u2264 0.75" = 2,
                              "\u2264 1.00" = 3)) +
    
    scale_size_identity() +
    scale_fill_manual(name = "Species", values = species_fills) +
    scale_color_manual(name = "Species", values = species_cols) +
    
    theme_graph() +
    theme(
      # legend.position = "inside",       # coordinates: x = left, y = top
      # legend.justification = c(0, 1),        # anchor top-left of legend box
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 10, face = "bold"),
      legend.key.size = unit(1, "lines")
    )
  
  return(graph_plot)
}

#########################
#     PLOT NETWORKS     #
#########################

# network including all species from the paper
allPorifera_graph_mtx <- get_graph_mtx_metacell(
  "05_SAMap_porifera/01_mapping_scores/AqueScilSlac_leiden3Clusters_100topCells_samapMappingTable.tsv",
  0.4)
allPorifera_plot <- plot_network_metacell(allPorifera_graph_mtx) +
  ggtitle("SAMap integration (threshold = 0.4)")
allPorifera_plot

ScilSlac_graph_mtx <- get_graph_mtx_metacell(
  "05_SAMap_porifera/01_mapping_scores/ScilSlac_leiden3Clusters_100topCells_samapMappingTable.tsv",
  0.4)
ScilSlac_plot <- plot_network_metacell(ScilSlac_graph_mtx) +
  ggtitle("SAMap integration (threshold = 0.4)")
ScilSlac_plot

AqueScil_graph_mtx <- get_graph_mtx_metacell(
  "05_SAMap_porifera/01_mapping_scores/AqueScil_leiden3Clusters_100topCells_samapMappingTable.tsv",
  0.4)
AqueScil_plot <- plot_network_metacell(AqueScil_graph_mtx) +
  ggtitle("SAMap integration (threshold = 0.4)")
AqueScil_plot

AqueSlac_graph_mtx <- get_graph_mtx_metacell(
  "05_SAMap_porifera/01_mapping_scores/AqueSlac_leiden3Clusters_100topCells_samapMappingTable.tsv",
  0.4)
AqueSlac_plot <- plot_network_metacell(AqueSlac_graph_mtx) +
  ggtitle("SAMap integration (threshold = 0.4)")
AqueSlac_plot

panel <- ggpubr::ggarrange(allPorifera_plot, ScilSlac_plot, AqueScil_plot, AqueSlac_plot,
                  nrow = 4)
panel

ggsave("05_SAMap_porifera/03_plots/samap_network_panel.pdf",
       panel, device = cairo_pdf,
       width = 8*1.5, height = 18*1.5, dpi = 300, bg = "white")
ggsave("05_SAMap_porifera/03_plots/samap_network_panel.png",
       panel, device = "png",
       width = 8*1.5, height = 18*1.5, dpi = 300, bg = "white")
