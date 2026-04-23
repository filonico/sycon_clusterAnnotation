#!/usr/bin/env Rscript

setwd("/lustre/alice3/data/evassvis/fn76/sycon/sycon_clusterAnnotation/sycon_cluster_analyses")

library(tidyverse)
library(ggsankeyfier)
library(ggpubr)
library(RColorBrewer)
library(purrr)


#############################
#     THEMESS FOR PLOTS     #
#############################

theme_for_plots <- theme(
  # aspect.ratio = 1,
  plot.background = element_rect(fill = "transparent", colour = NA), 
  panel.background = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_line(color = "grey90", lineend = "round"),
  panel.border = element_rect(colour = "black", linewidth = .6),
  legend.background = element_rect(fill = "transparent", colour = NA),
  legend.key = element_rect(fill = "transparent", colour = NA),
  legend.key.width = unit(.4, "cm"),
  legend.key.height = unit(.4, "cm"),
  legend.position = "right",
  legend.title = element_text(face = "bold", size = 10),
  legend.text = element_text(size = 10),
  axis.line = element_blank(),
  axis.ticks = element_line(colour = "black", linewidth = .4),
  axis.ticks.length = unit(0.10, "cm"),
  axis.text.x = element_text(color = "black",
                             margin = margin(t = 4, r = 0, b = 0, l = 0)),
  axis.text.y = element_text(color = "black",
                             margin = margin(t = 0, r = 4, b = 0, l = 0)),
  axis.title.y = element_text(angle = 90, size = 13,
                              margin = margin(t = 0, r = 10, b = 0, l = 0)),
  axis.title.x = element_text(angle = 0, size = 13,
                              margin = margin(t = 10, r = 0, b = 0, l = 0)),
  strip.text = element_text(color = "black", face = "bold", hjust = 0),
  strip.placement = "outside"
)


#####################
#     FUNCTIONS     #
#####################

# function to prepare the MappingScore SAMap matrix for plotting
tidyup_dataframe <- function(filename, threshold = 0.4) {
  
  df <- read.table(filename, header = TRUE, sep = "\t") %>%
    
    # create a long-format dataframe
    pivot_longer(-X, names_to = "target", values_to = "mapp_score") %>%
    rename("source" = "X") %>% 
    
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
    rename("connector" = "name",
           "node" = "value") %>%
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
  
  plot <- ggplot(data = df, aes(x = species, y = mapp_score,
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
    
    labs(x = "Species") +
    
    geom_text(data = . %>%
                filter(str_detect(connector, "from")),
              aes(label = node), cex = 3, hjust = 1, position = pos_text_left, stat = "sankeynode") +
    geom_text(data = . %>%
                filter(str_detect(connector, "to")),
              aes(label = node), cex = 3, hjust = 0, position = pos_text_right, stat = "sankeynode") +
    
    theme_bw(base_size = 12) +
    theme_for_plots +
    theme(panel.grid.major.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(face = "italic"))
  
  return(plot)
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


#########################################
#     SANKEYS SYCON VS DEMOSPONGIAE     #
#########################################

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
         species = case_when(species == "Aque" ~ "Amphimedon queenslandica",
                             species == "Scil" ~ "Sycon ciliatum",
                             TRUE ~ species),
         species = factor(species, levels = c("Sycon ciliatum", "Amphimedon queenslandica")),
         node = str_replace(node, "Scil_|Aque_", ""),
         node = str_replace_all(node, c("_" = " ",
                                        "Pinaco" = "Pinacocytes",
                                        "Archaeo" = "Archaeocytes",
                                        "Collagen" = "Collagen cells"))) %>%
  left_join(cluster_identity, by = join_by("node" == "cluster")) %>%
  mutate(annotation = str_replace_all(annotation, "_[0-9]+$", ""),
         annotation = str_replace_all(annotation, "_", " "))
# aquescil_mapping
aquescil_sankey <- aquescil_mapping %>%
  plot_sankey()
aquescil_sankey


scilslac_mapping <- tidyup_dataframe("05_SAMap_porifera/01_mapping_scores/ScilSlac_leiden3Clusters_100topCells_samapMappingTable.tsv") %>%
  mutate(species = case_when(species == "Slac" ~ "Spongilla lacustris",
                             species == "Scil" ~ "Sycon ciliatum",
                             TRUE ~ species),
         species = factor(species, levels = c("Sycon ciliatum", "Spongilla lacustris")),
         node = str_replace_all(node, c("Slac_10" = "Slac_10",
                                        "Slac_18" = "Slac_18",
                                        "Slac_15" = "Myopeptidocytes 1",
                                        "Slac_9" = "Myopeptidocytes 2",
                                        "Slac_8" = "Slac_8",
                                        "Slac_14" = "Metabolocytes 1",
                                        "Slac_32" = "Slac_32",
                                        "Slac_17" = "Pinacocytes",
                                        "Slac_24" = "Myopeptidocytes 3",
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

panel_sankeys <- ggpubr::ggarrange(aquescil_sankey, scilslac_sankey,
                                   common.legend = TRUE, legend = "right", labels = "AUTO")
panel_sankeys

ggsave("05_SAMap_porifera/03_plots/samap_sankey_panel.png",
       panel_sankeys, device = "png",
       width = 12, height = 6, dpi = 300, unit = "in", bg = "white")
ggsave("05_SAMap_porifera/03_plots/samap_sankey_panel.pdf",
       panel_sankeys, device = cairo_pdf,
       width = 12, height = 6, dpi = 300, unit = "in", bg = "white")

panel_sankeys <- ggpubr::ggarrange(aquescil_sankey + theme(axis.title.x = element_blank()),
                                   scilslac_sankey, nrow = 2,
                                   common.legend = TRUE, legend = "right")
panel_sankeys


#####################################
#     SANKEYS DEMOSPONGIAE ONLY     #
#####################################

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

# slacaque_mapping %>%
#   plot_sankey()


####################################
#     BOXPLOTS WITH CNIDARIANS     #
####################################

comparisons <- c("Porifera vs Cnidaria",
                 "Within Porifera",
                 "Within Cnidaria")

colors <- list(main  = setNames(c("#ffd166", "#ef476f", "#26547c"), comparisons),
               alpha = setNames(c("#f6e4ac", "#f98db3", "#8797a4"), comparisons))

# read in mapping scores
df_stitched_samap <- read.table("15_SAMap_cnidaria/01_mapping_scores/AqueHvulNvecScilSlacSpisXesp_leiden3Clusters_100topCells_samapMappingTable_leidenClusters.tsv",
                                header = TRUE, sep = "\t") %>%
  
  # create a long-format dataframe
  pivot_longer(-X, names_to = "target_ID", values_to = "mapp_score") %>%
  rename("source_ID" = "X") %>% 
  
  # filter out rows with mapping scores below the threshold
  filter(mapp_score >= 0.4) %>%
  
  # keep just one mapping of each pair   
  mutate(pair = paste(pmin(source_ID, target_ID), pmax(source_ID, target_ID), sep = "_")) %>%
  distinct(pair, .keep_all = TRUE) %>%
  select(source_ID, target_ID, mapp_score) %>%
  
  separate(source_ID, into = c("source_species", "source_cell"),
           extra = "merge", sep = "_", remove = FALSE) %>%
  separate(target_ID, into = c("target_species", "target_cell"),
           extra = "merge", sep = "_", remove = FALSE) %>%
  
  mutate(species_pair = paste0(source_species, "_", target_species),
         category_broad = case_when(species_pair %in% c("Aque_Scil", "Aque_Slac",
                                                        "Slac_Scil") ~ "Within Porifera",
                                    
                                    species_pair %in% c("Hvul_Nvec", "Hvul_Xesp",
                                                        "Hvul_Spis", "Nvec_Xesp",
                                                        "Nvec_Spis", "Xesp_Hvul",
                                                        "Xesp_Nvec", "Xesp_Spis",
                                                        "Spis_Hvul", "Spis_Nvec",
                                                        "Spis_Xesp") ~ "Within Cnidaria",
                                    
                                    TRUE ~ "Porifera vs Cnidaria"), 
         category_broad = as_factor(category_broad),
         category_broad = relevel(category_broad, ref = "Within Porifera")) %>%
  arrange(mapp_score) %>%
  group_by(category_broad) %>%
  mutate(rank = row_number(),
         survival = 1 - rank / n(),
         # survival = case_when(survival == 0 ~ 2e-16,
         #                      TRUE ~ survival)
         high_score = mapp_score > 0.9)


wilcox_comparisons <- combn(levels(df_stitched_samap$category_broad), 2, simplify = FALSE)
pairwise.wilcox.test(df_stitched_samap$mapp_score, df_stitched_samap$category_broad, p.adj = "bonferroni")

cnidaria_scores <- subset(df_stitched_samap, category_broad == "Within Cnidaria")$mapp_score
porifera_scores <- subset(df_stitched_samap, category_broad == "Within Porifera")$mapp_score
cnidariaVSporifera_scores <- subset(df_stitched_samap, category_broad == "Porifera vs Cnidaria")$mapp_score

ks.test(porifera_scores, cnidaria_scores, alternative = "greater")
ks.test(porifera_scores, cnidariaVSporifera_scores, alternative = "greater")
effsize::cliff.delta(cnidaria_scores, porifera_scores)

boxplot <- df_stitched_samap %>%
  ggplot(aes(x = mapp_score, y = category_broad, col = category_broad, fill = category_broad)) +
  geom_jitter(size = 2, height = 0.4) +
  
  scale_color_manual(values = colors$alpha) +

  geom_boxplot(width = 0.4, col = "black", linewidth = 0.7,
               outliers = FALSE, staplewidth = 0.5) +
  
  # stat_compare_means(comparisons = wilcox_comparisons, method = "wilcox.test",
  #                    symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
  #                                       symbols = c("****", "***", "**", "*", "ns"))) +
  
  scale_fill_manual(values = colors$main) +
  scale_y_discrete(limits = rev,
                   labels = c(paste0("Within\nCnidaria\n(n = ", length(cnidaria_scores), ")"),
                              paste0("Porifera vs\nCnidaria\n(n = ", length(cnidariaVSporifera_scores), ")"),
                              paste0("Within\nPorifera\n(n = ", length(porifera_scores), ")"))) +

  labs(x = "Mapping scores", y = "Comparison") +

  theme_bw(base_size = 12) +
  theme_for_plots +
  theme(legend.position = "none",
        axis.title.y = element_blank())

boxplot

survival_plot <- df_stitched_samap %>%
  ggplot(aes(x = mapp_score, y = survival)) +
  
  geom_line(aes(color = category_broad), linewidth = 1) +
  scale_color_manual(values = colors$main) +
  
  # geom_text(x = 0.55, y = 0.9, label = "KS test, p = 0.002") +
  
  labs(x = "Mapping score", y = "P(score ≥ x)")  +
  
  guides(color = guide_legend(keyheight = 0.8,
                              default.unit = "cm")) +
  
  theme_bw(base_size = 12) +
  theme_for_plots +
  theme(legend.position = "inside",
        legend.position.inside = c(0.72, 0.82),
        legend.background = element_rect(fill = "white",
                                         colour = "black", linewidth = .4),
        legend.title = element_blank(),
        legend.margin = margin(0, 6, 0, 3, "mm"))
survival_plot

panel_statistics <- ggarrange(boxplot + theme(axis.title.x = element_blank()),
                              survival_plot,
                              nrow = 2) +
  theme(plot.margin = margin(0, 0, 0, 10, "mm"))

final_panel <- ggarrange(panel_sankeys, panel_statistics,
                         labels = "AUTO", font.label = list(size = 20),
                         align = "hv", widths = c(1, 0.7))
final_panel

ggsave("05_SAMap_porifera/03_plots/final_panel_fig3.png",
       final_panel, device = "png",
       width = 14/1.1, height = 10/1.1, dpi = 300, unit = "in", bg = "white")
ggsave("05_SAMap_porifera/03_plots/final_panel_fig3.pdf",
       final_panel, device = cairo_pdf,
       width = 14/1.1, height = 10/1.1, dpi = 300, unit = "in", bg = "white")
