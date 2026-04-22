#!/usr/bin/env Rscript

setwd("/lustre/alice3/data/evassvis/fn76/sycon/sycon_clusterAnnotation/sycon_cluster_analyses")

library(tidyverse)
library(ggsankeyfier)
library(ggpubr)
library(RColorBrewer)
library(purrr)


#####################
#     FUNCTIONS     #
#####################

# function to prepare the MappingScore SAMap matrix for plotting
tidyup_dataframe <- function(filename, threshold = 0.4) {
  
  df <- read.table(filename, header = TRUE, sep = "\t") %>%
    
    # create a long-format dataframe
    pivot_longer(-X, names_to = "target", values_to = "mapp_score") %>%
    rename("X" = "source") %>% 
    
    # filter out rows with mapping scores below the threshold
    filter(mapp_score >= threshold) %>%
    
    # keep just one mapping of each pair   
    mutate(pair = paste(pmin(source, target), pmax(source, target), sep = "_")) %>%
    distinct(pair, .keep_all = TRUE) %>%
    select(source, target, mapp_score) %>%
    
    mutate(edge_id = paste0(source, "_", target)) %>%
    pivot_longer(-c(mapp_score, edge_id)) %>%
    mutate(edge_id = as.integer(factor(edge_id)),
           name = str_replace(name, "source", "from"),
           name = str_replace(name, "target", "to")) %>%
    rename("name" = "connector",
           "value" = "node") %>%
    separate(node, into = c("species", "cell"),
             extra = "merge", sep = "_", remove = FALSE) %>%
    select(-cell) %>%
    
    # create a column with discrete categories for mapping scores
    mutate(category = factor(cut(mapp_score, breaks = value_intervals,
                                 labels = as.factor(value_interval_labels),
                                 right = FALSE),
                             levels = c("\u2264 1.00", "\u2264 0.75", "\u2264 0.50", "\u2264 0.25"))) 
  
  return(df)
}


# function to produce a sankey plot out of the tidied-up dataframe
plot_sankey <- function(df) {
  
  ggplot(data = df, aes(x = species, y = mapp_score,
             group = node, connector = connector, edge_id = edge_id)) +
    
    geom_point(aes(fill = category), x = Inf, y = Inf,
               col = NA, alpha = 0.5,
               size = 0, shape = 22) +
    ggsankeyfier::geom_sankeyedge(aes(fill = category), alpha = 0.5,
                                  position = pos, show.legend = FALSE) +
    
    scale_fill_grey(name = "Mapping\nscore",
                    start = 0.1, end = 0.8) +
    guides(fill = guide_legend(override.aes = list(size = 3))) +
    
    ggnewscale::new_scale_fill() +
    ggnewscale::new_scale_fill() +
    
    geom_point(aes(fill = annotation), col = "grey20", x = Inf, y = Inf, size = 0, shape = 22) +
    ggsankeyfier::geom_sankeynode(aes(fill = annotation), col = "grey20",
                                  linewidth = 0.4, position = pos, show.legend = FALSE) +
    
    scale_fill_manual(values = c("Accessory cells" = "#B49440",
                                 "Archaeocyte-like\nstem cells" = "#B65CBF",
                                 "Choanocytes" = "#67A64E",
                                 "Metabolocyte- and\npinacocyte-like\nearly embryos" = "#747CC9",
                                 "Myopeptidocyte-like\nchoanocytes" = "#CA5F3E",
                                 "Oocytes/\nearly embryos" = "#48B1A7",
                                 "Sclerocyte-like\npinacocytes" = "#C8577B",
                                 "Uncertain" = "grey60",
                                 "Unknown" = "grey85",
                                 "Other species" = "white"),
                      na.value = "white",
                      breaks = ~ .x[!is.na(.x)],
                      name = "Cluster\nannotation") +
    guides(fill = guide_legend(override.aes = list(size = 3),
                               keyheight = 1.1,
                               default.unit = "cm")) +
    
    geom_text(data = . %>%
                filter(str_detect(connector, "from")),
              aes(label = node), cex = 3, hjust = 1, position = pos_text_left, stat = "sankeynode") +
    geom_text(data = . %>%
                filter(str_detect(connector, "to")),
              aes(label = node), cex = 3, hjust = 0, position = pos_text_right, stat = "sankeynode") +
    
    theme_bw(base_size = 12) +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_text(face = "italic", margin = margin(8, 0, 0, 0)),
          axis.text.y = element_blank(),
          legend.title = element_text(face = "bold", size = 10),
          legend.text = element_text(size = 10))
}



get_mappingScore_distribution <- function(filename, experiment = "leiden3Clusters") {
  
  dataframe <- read.table(filename, header = TRUE, sep = "\t", row.names = 1)

  dataframe[upper.tri(dataframe)] <- NA

  dataframe <- dataframe %>%
    rownames_to_column(var = 'source') %>%
    pivot_longer(-source, names_to = "target", values_to = "value") %>%
    drop_na() %>%
    separate(source, into = c("source_species", "source_group"), sep = "_", extra = "merge") %>%
    separate(target, into = c("target_species", "target_group"), sep = "_", extra = "merge") %>%
    filter(source_species != target_species) %>%
    mutate(comparison = paste0(source_species, "_", target_species),
           experiment = experiment) %>%
    select(c(comparison, experiment, value))
  
  return(dataframe)
}


############################
#     GENERATE SANKEYS     #
############################

# define grey ranges for sankey edges
value_intervals <- c(0, 0.25, 0.5, 0.75, 1)
value_interval_labels <- c("\u2264 0.25", "\u2264 0.50", "\u2264 0.75", "\u2264 1.00")

# define positions of sankey elements
pos <- position_sankey(v_space = "auto", order = "ascending")
pos_text_left <- position_sankey(v_space = "auto", order = "ascending", nudge_x = -0.1)
pos_text_right <- position_sankey(v_space = "auto", order = "ascending", nudge_x = 0.1)

# define sycon cluster annotation
cluster_identity <- tibble(cluster = seq(0, 32),
                           annotation = c("Uncertain_0", "Choanocytes_1", "Unknown_2", "Unknown_3",
                                          "Unknown_4", "Uncertain_5", "Unknown_6", "Choanocytes_7", "Unknown_8",
                                          "Unknown_9", "Uncertain_10", "Choanocytes_11", "Unknown_12", "Unknown_13",
                                          "Unknown_14", "Metabolocyte-_and\npinacocyte-like\nearly_embryos",
                                          "Choanocytes_16", "Oocytes/\nearly_embryos", "Unknown_18", "Oocytes/\nearly_embryos",
                                          "Unknown_20", "Myopeptidocyte-like\nchoanocytes", "Unknown_22", "Choanocytes_23",
                                          "Oocytes/\nearly_embryos", "Sclerocyte-like\npinacocytes", "Uncertain_26",
                                          "Oocytes/\nearly_embryos", "Unknown_28", "Unknown_29", "Accessory_cells",
                                          "Unknown_31", "Archaeocyte-like\nstem_cells")) %>%
  mutate(cluster = as.character(cluster))
cluster_identity

aquescil_mapping <- tidyup_dataframe("05_SAMap_porifera/01_mapping_scores/AqueScil_leiden3Clusters_100topCells_samapMappingTable.tsv") %>%
  mutate(connector = case_when(connector == "from" ~ "to",
                               connector == "to" ~ "from",
                               TRUE ~ connector),
         species = case_when(species == "Aque" ~ "Amphimedon\nqueenslandica",
                             species == "Scil" ~ "Sycon\nciliatum",
                             TRUE ~ species),
         species = factor(species, levels = c("Sycon\nciliatum", "Amphimedon\nqueenslandica")),
         node = str_replace(node, "Scil_|Aque_", ""),
         node = str_replace_all(node, c("_" = " ",
                                        "Pinaco" = "Pinacocytes",
                                        "Archaeo" = "Archaeocytes",
                                        "Collagen" = "Collagen cells"))) %>%
  left_join(cluster_identity, by = join_by("node" == "cluster")) %>%
  mutate(annotation = str_replace_all(annotation, "_[0-9]+$", ""),
         annotation = str_replace_all(annotation, "_", " "))
aquescil_mapping
aquescil_sankey <- aquescil_mapping %>%
  plot_sankey()
aquescil_sankey


scilslac_mapping <- tidyup_dataframe("05_SAMap_porifera/01_mapping_scores/ScilSlac_leiden3Clusters_100topCells_samapMappingTable.tsv") %>%
  mutate(species = case_when(species == "Slac" ~ "Spongilla\nlacustris",
                             species == "Scil" ~ "Sycon\nciliatum",
                             TRUE ~ species),
         species = factor(species, levels = c("Sycon\nciliatum", "Spongilla\nlacustris")),
         node = str_replace_all(node, c("Slac_10" = "Slac_10",
                                        "Slac_18" = "Slac_18",
                                        "Slac_15" = "Myopeptido-\ncytes 1",
                                        "Slac_9" = "Myopeptido-\ncytes 2",
                                        "Slac_8" = "Slac_8",
                                        "Slac_14" = "Metabolocytes 1",
                                        "Slac_32" = "Slac_32",
                                        "Slac_17" = "Pinacocytes",
                                        "Slac_24" = "Myopeptido-\ncytes 3",
                                        "Slac_13" = "Slac_13",
                                        "Slac_0" = "Archaeocytes 1",
                                        "Slac_22" = "Slac_22",
                                        "Slac_25" = "Choanocytes",
                                        "Slac_4" = "Archaeocytes 5")),
         node = str_replace(node, "Scil_", "")) %>%
  left_join(cluster_identity, by = join_by("node" == "cluster")) %>%
  mutate(annotation = str_replace(annotation, "_[0-9]+$", ""),
         annotation = str_replace_all(annotation, "_", " "),
         node = str_replace(node, "Slac_", ""))
scilslac_mapping
scilslac_sankey <- scilslac_mapping %>%
  plot_sankey()
scilslac_sankey

panel <- ggpubr::ggarrange(aquescil_sankey, scilslac_sankey,
                           common.legend = TRUE, legend = "right", labels = "AUTO")

ggsave("05_SAMap_porifera/03_plots/samap_sankey_panel.png",
       panel, device = "png",
       width = 12, height = 6, dpi = 300, unit = "in", bg = "white")
ggsave("05_SAMap_porifera/03_plots/samap_sankey_panel.pdf",
       panel, device = cairo_pdf,
       width = 12, height = 6, dpi = 300, unit = "in", bg = "white")


slacaque_mapping <- tidyup_dataframe("05_SAMap_porifera/01_mapping_scores/AqueSlac_leiden3Clusters_100topCells_samapMappingTable.tsv") %>%
  mutate(species = case_when(species == "Slac" ~ "S. lacustris",
                             species == "Aque" ~ "A. queenslandica",
                             TRUE ~ species),
         species = factor(species, levels = c("A. queenslandica", "S. lacustris")),
         node = str_replace_all(node, c("Slac_10" = "10",
                                        "Slac_18" = "18",
                                        "Slac_15" = "Myopeptidocytes 1",
                                        "Slac_9" = "Myopeptidocytes 2",
                                        "Slac_8" = "8",
                                        "Slac_14" = "Metabolocytes 1",
                                        "Slac_32" = "32",
                                        "Slac_17" = "Pinacocytes",
                                        "Slac_24" = "Myopeptidocytes 3",
                                        "Slac_13" = "13",
                                        "Slac_0" = "Archaeocytes 1",
                                        "Slac_22" = "22",
                                        "Slac_25" = "Choanocytes/-blasts",
                                        "Slac_4" = "Archaeocytes 5")),
         node = str_replace_all(node, c("Aque_Pinaco" = "Pinacocytes",
                                        "Aque_Archaeo" = "Archaeocytes",
                                        "Aque_Collagen" = "Collagen cells",
                                        "Aque_Choano_to_pinaco" = "Choano- to pinacocytes",
                                        "Aque_Choanocytes" = "Choanocytes",
                                        "Aque_Sperm" = "Sperm",
                                        "Aque_Aspcizin" = "Aspcizin",
                                        "Aque_Bactericidial" = "Bactericidial cells",
                                        "Aque_Unk" = "Unknown",
                                        "_" = " ")),
         node = str_replace(node, "Aque |Slac ", ""))

slacaque_mapping %>%
  plot_sankey()

aquescil_mapping %>%
  add_row(scilslac_mapping %>%
            mutate(edge_id = edge_id + nrow(aquescil_mapping)/2)) %>%
  plot_sankey()


ScilSlac_sankey <- fromTable_toSankey(
  "05_SAMap_porifera/01_mapping_scores/ScilSlac_leiden3Clusters_100topCells_samapMappingTable.tsv",
  "Scil", "Slac", "Sycon", "Spongilla", 0.4
)
ScilSlac_sankey

AqueScil_sankey <- fromTable_toSankey(
  "05_SAMap_porifera/01_mapping_scores/AqueScil_leiden3Clusters_100topCells_samapMappingTable.tsv",
  "Scil", "Aque", "Sycon", "Amphimedon", 0.4
)
AqueScil_sankey

AqueSlac_sankey <- fromTable_toSankey(
  "05_SAMap_porifera/01_mapping_scores/AqueSlac_leiden3Clusters_100topCells_samapMappingTable.tsv",
  "Slac", "Aque", "Spongilla", "Amphimedon", 0.4
)
AqueSlac_sankey

panel <- ggpubr::ggarrange(ScilSlac_sankey, AqueScil_sankey, AqueSlac_sankey,
                  nrow = 3, common.legend = TRUE, legend = "right")
panel

# save the sankey
ggsave("05_SAMap_porifera/03_plots/samap_sankey_panel.pdf",
     plot = panel, device = cairo_pdf,
     dpi = 300, height = 14, width = 8, units = ("in"), bg = 'white')
ggsave("05_SAMap_porifera/03_plots/samap_sankey_panel.png",
       plot = panel, device = "png",
       dpi = 300, height = 14, width = 8, units = ("in"), bg = 'white')


ScilBlobSlac_sankey <- fromTable_toSankey(
  "13_recluster_blob/04_SAMap/ScSlac_leiden3Clusters_0topCells_samapMappingTable.tsv",
  "Scil", "Slac", "Sycon", "Spongilla", 0.2
)
