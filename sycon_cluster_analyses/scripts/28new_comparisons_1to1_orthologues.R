#!/usr/bin/env Rscript

setwd("/data/evassvis/fn76/sycon/sycon_clusterAnnotation/sycon_cluster_analyses/")

library(tidyverse)
library(Seurat)
library(SeuratExtend)
library(VennDiagram)
library(patchwork)

###############################
#     DEFINE PLOT THEMES      #
###############################

theme_for_plots <- theme(
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
  axis.text.x = element_text(color = "black", margin = margin(t = 4, r = 0, b = 0, l = 0)),
  axis.text.y = element_text(color = "black", margin = margin(t = 0, r = 4, b = 0, l = 0)),
  axis.title.y = element_text(angle = 90, size = 13, margin = margin(t = 0, r = 10, b = 0, l = 0)),
  axis.title.x = element_text(angle = 0, size = 13, margin = margin(t = 10, r = 0, b = 0, l = 0)),
  strip.text = element_text(color = "black", face = "bold", hjust = 0),
  strip.placement = "outside",
  strip.background = element_blank(),
  strip.clip = "off"
)

umap_arrows <- list(
  annotation_custom(grob = grid::segmentsGrob(x0 = unit(0, "mm"), x1 = unit(12, "mm"),
                                              y0 = unit(0, "mm"), y1 = unit(0, "mm"),
                                              arrow = arrow(length = unit(2.5, "mm"), ends = "last", type = "open"),
                                              gp = grid::gpar(col = "black", fill = "black", lwd = 1))),
  annotation_custom(grob = grid::segmentsGrob(x0 = unit(0, "mm"), x1 = unit(0, "mm"),
                                              y0 = unit(0, "mm"), y1 = unit(12, "mm"),
                                              arrow = arrow(length = unit(2.5, "mm"), ends = "last", type = "open"),
                                              gp = grid::gpar(col = "black", fill = "black", lwd = 1))),
  annotation_custom(grob = grid::textGrob(label = "UMAP 1", x = unit(0, "mm"), y = unit(0, "mm") - unit(2.5, "mm"),
                                          just = c(0, 1), gp = grid::gpar(fontsize = 10))),
  annotation_custom(grob = grid::textGrob(label = "UMAP 2", x = unit(0, "mm") - unit(2.5, "mm"), y = unit(0, "mm"),
                                          just = c(0, 0), rot = 90, gp = grid::gpar(fontsize = 10))),
  coord_cartesian(clip = "off")
)

theme_for_UMAPS <- theme(
  plot.background = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  panel.grid = element_blank(),
  legend.text = element_text(size = 10),
  legend.title = element_text(size = 10, face = "bold"),
  plot.title = element_text(size = 13, hjust = 0.5, vjust = 1.75, face = "bold"),
  axis.line = element_blank(),
  axis.ticks = element_blank(),
  axis.text = element_blank(),
  axis.title = element_blank()
)


#####################
#     LOAD DATA     #
#####################

PCA_1to1 <- read.table("17_comparisons_1to1_orthologues/01_orthologue_tables/1_to_1_Sycon_Spongilla.tsv",
                       header = TRUE, sep = "\t", na.strings = c("", ".")) %>%
  select(-c(required_nodes_present, required_single_nodes_present, collected_specs_present)) %>%
  as_tibble()

MCA_1to1 <- read.table("17_comparisons_1to1_orthologues/01_orthologue_tables/1_to_1_Sycon_Spongilla_present_on_the_metazoan_stem.tsv",
                       header = TRUE, sep = "\t", na.strings = c("", ".")) %>%
  select(-c(required_nodes_present, required_single_nodes_present, collected_specs_present)) %>%
  as_tibble()

PCA_1toMany <- read.table("17_comparisons_1to1_orthologues/01_orthologue_tables/1_to_many_Sycon_Spongilla.tsv",
                          header = TRUE, sep = "\t", na.strings = c("", ".")) %>%
  select(-c(required_nodes_present, required_single_nodes_present, collected_specs_present)) %>%
  as_tibble()

MCA_1toMany <- read.table("17_comparisons_1to1_orthologues/01_orthologue_tables/1_to_many_Sycon_Spongilla_present_on_the_metazoan_stem.tsv",
                          header = TRUE, sep = "\t", na.strings = c("", ".")) %>%
  select(-c(required_nodes_present, required_single_nodes_present, collected_specs_present)) %>%
  as_tibble()

# check overlaps
venn.diagram(list("PCA_1to1" = PCA_1to1$orthogroup_id,
                  "MCA_1to1" = MCA_1to1$orthogroup_id,
                  "PCA_1toMany" = PCA_1toMany$orthogroup_id,
                  "MCA_1toMany" = MCA_1toMany$orthogroup_id),
             filename = NULL, fill = c("red", "green", "blue", "yellow"))


# slac conversion table
slac_gene_names_conversion <- read.table("../spongilla_remapping/00_input/slac_genome/genes_to_transcript.tsv",
                                         sep = "\t", na.strings = "N/A") %>%
  left_join(read.table("../spongilla_remapping/00_input/slac_genome/transcripts_to_cdss.tsv",
                       sep = "\t", na.strings = "N/A"),
            by = join_by("V1" == "V2")) %>%
  as_tibble() %>%
  rename("transcript_id" = "V1", "gene_id" = "V2", "protein_id" = "V1.y") %>%
  arrange(gene_id)

# integrated samap
ScilSlac_samap <- schard::h5ad2seurat("05NEW_SAMap_porifera/ScilSlac_leiden3Clusters_samap.h5ad")
ScilSlac_samap[[]]$Scil_leiden_clusters <- as.factor(ScilSlac_samap[[]]$Scil_leiden_clusters)
ScilSlac_samap[[]]$Slac_leiden_clusters <- as.factor(ScilSlac_samap[[]]$Slac_leiden_clusters)
ScilSlac_samap[[]]$leiden_clusters      <- as.factor(ScilSlac_samap[[]]$leiden_clusters)


###############################################
#     COUNT HOW MANY GENES PER ORTHOGROUP     #
###############################################

# create a single df with all counts per orthogroup per species
orthogroup_full_set <- bind_rows(PCA_1to1, PCA_1toMany, MCA_1to1, MCA_1toMany) %>%
  # keep only one line per OG
  distinct(orthogroup_id, .keep_all = TRUE) %>%
  pivot_longer(c(LDemCA__Spongilla_lacustris, LCalcCA__Sycon_ciliatum),
               names_to = "species", values_to = "gene_ID") %>%
  # fix species anem
  mutate(species = case_when(str_detect(species, "Spongilla") ~ "Slac",
                             str_detect(species, "Sycon") ~ "Scil",
                             TRUE ~ species)) %>%
  # split multi-copy OG into multiple lines
  separate_rows(gene_ID, sep = ";") %>%
  # adjust slac gene names
  mutate(gene_ID = str_remove(gene_ID, "\\.[0-9]+$")) %>%
  left_join(slac_gene_names_conversion, by = join_by("gene_ID" == "protein_id")) %>%
  # add gene names for slac
  mutate(gene_id = if_else(is.na(gene_id), gene_ID, gene_id)) %>%
  select(-c(transcript_id, gene_ID)) %>%
  # remove duplicated genes (coming from slac isoforms)
  distinct()

# count genes per OG
orthogroup_full_set_with_geneCounts <- orthogroup_full_set %>%
  group_by(orthogroup_id, species) %>%
  mutate(n_genes = n()) %>%
  ungroup()


#######################################
#     FIND MARKERS PER SPECIES        #
#######################################

# skip this block if rerunning and marker genes can be input from tables

# subset sycon and SCTransform
scil_samap <- ScilSlac_samap %>% subset(species == "Scil") %>%
  SCTransform()

# find markers on sycon specific leiden clusters
scil_samap_allMarkers <- scil_samap %>%
  FindAllMarkers(group.by = "Scil_leiden_clusters") %>%
  mutate(gene = str_replace(gene, "Scil-", ""))

write.table(scil_samap_allMarkers,
            "17_comparisons_1to1_orthologues/scil_samap_allMarkers_leidenClusters.tsv",
            col.names = TRUE, sep = "\t", quote = FALSE, row.names = FALSE)

# subset spongilla and SCTransform
slac_samap <- ScilSlac_samap %>% subset(species == "Slac") %>%
  SCTransform()

# find markers on spongilla specific leiden clusters
slac_samap_allMarkers <- slac_samap %>%
  FindAllMarkers(group.by = "Slac_leiden_clusters") %>%
  mutate(gene = str_replace(gene, "Slac-", "")) %>%
  # add gene_id conversion
  left_join(slac_gene_names_conversion, by = join_by("gene" == "gene_id"),
            relationship = "many-to-many") %>%
  select(-c(transcript_id, protein_id))

write.table(slac_samap_allMarkers,
            "17_comparisons_1to1_orthologues/slac_samap_allMarkers_leidenClusters.tsv",
            col.names = TRUE, sep = "\t", quote = FALSE, row.names = FALSE)


#####################################################################
#     COMPUTE NUMBER OF CLUSTERS WHERE EACH OG IS EXPRESSED         #
#####################################################################

scil_samap_allMarkers <- read.table("17_comparisons_1to1_orthologues/scil_samap_allMarkers_leidenClusters.tsv",
                                    header = TRUE, sep = "\t") %>% as_tibble()
slac_samap_allMarkers <- read.table("17_comparisons_1to1_orthologues/slac_samap_allMarkers_leidenClusters.tsv",
                                    header = TRUE, sep = "\t") %>% as_tibble()

# create a single dataframe with all markers per cluster
allMarkers_per_cluster <- bind_rows(scil_samap_allMarkers %>%
                                      mutate(cluster = paste0("Scil_", cluster)),
                                    slac_samap_allMarkers %>%
                                      mutate(cluster = paste0("Slac_", cluster)) %>%
                                      distinct()) %>%
  filter(p_val_adj < 0.05, avg_log2FC > 1) %>%
  select(cluster, gene)
allMarkers_per_cluster

# add cluster-expression info at the gene level
orthogroup_full_set_with_clusterCount <- orthogroup_full_set %>%
  left_join(allMarkers_per_cluster, by = join_by("gene_id" == "gene"),
            relationship = "many-to-many") %>%
  group_by(orthogroup_id, species) %>%
  mutate(n_clusters = n_distinct(cluster, na.rm = TRUE)) %>%
  ungroup()

# compute og_type classification
og_type_table <- orthogroup_full_set_with_geneCounts %>%
  distinct(orthogroup_id, species, n_genes) %>%
  pivot_wider(names_from = species, values_from = n_genes) %>%
  mutate(og_type = case_when(Scil == 1 & Slac == 1 ~ "1 to 1",
                             Scil > 1  & Slac == 1 ~ "Single copy\nin Spongilla",
                             Scil == 1 & Slac > 1 ~ "Single copy\nin Sycon",
                             TRUE ~ "Many to many")) %>%
  select(orthogroup_id, og_type)

# get the total number of per-species leiden clusters
n_clusters_total <- tibble(species = c("Scil", "Slac"),
                           n_clusters_total = c(nlevels(droplevels(ScilSlac_samap[[]]$Scil_leiden_clusters)),
                                                nlevels(droplevels(ScilSlac_samap[[]]$Slac_leiden_clusters))))
n_clusters_total

# combine everything
orthogroup_full_set_with_geneCounts_clusterCounts <- orthogroup_full_set_with_geneCounts %>%
  left_join(orthogroup_full_set_with_clusterCount %>%
              select(orthogroup_id, species, cluster, n_clusters) %>%
              distinct(),
            by = c("orthogroup_id", "species"),
            relationship = "many-to-many") %>%
  left_join(og_type_table, by = "orthogroup_id") %>%
  left_join(n_clusters_total, by = "species") %>%
  mutate(n_cluster_prop = n_clusters / n_clusters_total) %>%
  group_by(orthogroup_id) %>%
  # remove all OGs where at least one gene is expressed in any cluster
  filter(!any(n_clusters == 0)) %>%
  ungroup() %>%
  mutate(species_label = case_when(species == "Scil" ~ "S. ciliatum",
                                   species == "Slac" ~ "S. lacustris"),
         og_type = factor(og_type, levels = c("1 to 1", "Single copy\nin Sycon",
                                              "Single copy\nin Spongilla", "Many to many")))

orthogroup_full_set_with_geneCounts_clusterCounts


#############################################################
#     GENE-LEVEL BREADTH METRICS (participation/redundancy) #
#############################################################

gene_breadth <- orthogroup_full_set_with_clusterCount %>%
  group_by(orthogroup_id, species, gene_id) %>%
  summarise(genewise_n_clusters = n_distinct(cluster, na.rm = TRUE),
            .groups = "drop")
gene_breadth

og_breadth_metrics <- gene_breadth %>%
  group_by(orthogroup_id, species) %>%
  summarise(n_genes = n(),
            n_genes_asMarkers = sum(genewise_n_clusters > 0),
            n_marker_events = sum(genewise_n_clusters),
            .groups = "drop") %>%
  left_join(orthogroup_full_set_with_geneCounts_clusterCounts %>%
              distinct(orthogroup_id, species, n_clusters, n_cluster_prop, og_type),
            by = c("orthogroup_id", "species")) %>%
  mutate(
    # fraction of an OG's genes that are a marker of *something*
    # (doesn't distinguish overlap vs. distinct clusters -- see redundancy)
    participation = n_genes_asMarkers / n_genes,
    # total marker-events / distinct clusters marked; ==1 no overlap between
    # paralogs, >1 some clusters marked by multiple paralogs, 0 no marker anywhere
    redundancy    = n_marker_events / pmax(n_clusters, 1)
  )

og_breadth_metrics


#####################################
#     PAIRED SPECIES COMPARISON     #
#####################################

breadth_wide <- og_breadth_metrics %>%
  select(orthogroup_id, species, og_type, n_cluster_prop) %>%
  pivot_wider(names_from = species, values_from = n_cluster_prop) %>%
  drop_na(Scil, Slac, og_type)

# wide format guarantees pairing for wilcoxon paired test
breadth_long <- breadth_wide %>%
  pivot_longer(c(Scil, Slac), names_to = "species", values_to = "n_cluster_prop") %>%
  mutate(species_label = if_else(species == "Scil", "S. ciliatum", "S. lacustris"))

# paired Wilcoxon signed-rank test per og_type
paired_tests_breadth <- breadth_wide %>%
  group_by(og_type) %>%
  summarise(n_pairs = n(),
            p_value = wilcox.test(Scil, Slac, paired = TRUE)$p.value,
            .groups = "drop") %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"),
         label = case_when(p_adj < 0.001 ~ "***",
                           p_adj < 0.01 ~ "**",
                           p_adj < 0.05 ~ "*",
                           TRUE ~ "ns"))
paired_tests_breadth

count_og_type <- breadth_long %>%
  count(og_type) %>%
  mutate(n = n/2) %>%
  deframe()

boxplots_per_og_type <- breadth_long %>%
  ggplot(aes(x = species_label, y = n_cluster_prop)) +
  geom_line(aes(group = orthogroup_id), alpha = 0.4, colour = "grey70") +
  geom_jitter(aes(col = species_label), width = 0.15, height = 0.01) +
  geom_boxplot(outliers = FALSE, fill = alpha("white", 0),
               linewidth = 0.7, width = 0.4) +
  
  geom_text(data = paired_tests_breadth,
            aes(x = 1.5, y = 0.32, label = label),
            inherit.aes = FALSE, size = 4) +
  # annotate(geom = "segment", x = 1, xend = 2, y = 0.312, linewidth = 0.5) +
  
  scale_y_continuous(limits = c(0, 0.33)) +
  scale_color_manual(values = c("S. ciliatum" = "#ff9ebb",
                                "S. lacustris" = "#b9375e"),
                     guide = "none") +
  
  labs(y = "Proportion of clusters") +
  
  facet_wrap(~og_type, ncol = 4,
             labeller = labeller(og_type = as_labeller(
               c("1 to 1" =  paste0("1 to 1\n(n = ", count_og_type[["1 to 1"]], ")"),
                 "Single copy\nin Sycon" =  paste0("Single copy in Sycon\n(n = ", count_og_type[["Single copy\nin Sycon"]], ")"),
                 "Single copy\nin Spongilla" =  paste0("Single copy in Spongilla\n(n = ", count_og_type[["Single copy\nin Spongilla"]], ")"),
                 "Many to many" =  paste0("Many to many\n(n = ", count_og_type[["Many to many"]], ")"))))) +
  theme_bw(base_size = 12) +
  theme_for_plots +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, face = "italic", hjust = 1))
boxplots_per_og_type

# # ridge distributions (unchanged in spirit, species_label used for consistency)
# ridges_per_species <- orthogroup_full_set_with_geneCounts_clusterCounts %>%
#   ggplot(aes(x = n_cluster_prop, y = species_label, fill = species_label, col = species_label)) +
#   ggridges::geom_density_ridges(bandwidth = 0.02, alpha = 0.8, linewidth = 0.6) +
#   scale_color_manual(values = c("S. ciliatum" = "#ff9ebb", "S. lacustris" = "#b9375e"), guide = "none") +
#   scale_fill_manual(values = c("S. ciliatum" = "#ff9ebb", "S. lacustris" = "#b9375e"), guide = "none") +
#   scale_x_continuous(breaks = seq(0, 0.3, 0.1), limits = c(0, 0.35)) +
#   coord_flip() +
#   facet_wrap(~species_label) +
#   theme_bw(base_size = 12) +
#   theme_for_plots +
#   theme(panel.border = element_blank(),
#         axis.title.x = element_blank(), axis.title.y = element_blank(),
#         axis.text.x = element_blank(), axis.ticks.x = element_blank())
# ridges_per_species
# 
# panel <- (boxplots_per_og_type + ridges_per_species) +
#   plot_layout(axes = "collect", widths = c(1, 0.5))
# panel
# 
ggsave("17_comparisons_1to1_orthologues/boxplot_orthogroup_categories.png",
       boxplots_per_og_type, device = "png",
       height = 4, width = 8, dpi = 300, unit = "in", bg = "white")
ggsave("17_comparisons_1to1_orthologues/boxplot_orthogroup_categories.pdf",
       boxplots_per_og_type, device = cairo_pdf,
       height = 4, width = 8, dpi = 300, unit = "in", bg = "white")


###################################################################
#     PAIRED SPECIES COMPARISON MANY-TO-MANY ORTHOGROUPS ONLY     #
###################################################################

# Only many-to-many OGs have a genuinely comparable, non-degenerate
# participation/redundancy value on BOTH sides -- 1-to-1 is trivially 1 gene
# (participation always 1, redundancy undefined), and the two single-copy
# categories only have a non-trivial value in the species that duplicated.
# See mtm_metrics / single_copy_metrics below for how those are handled.

mtm_wide_participation <- og_breadth_metrics %>% 
  filter(og_type != "1 to 1") %>%
  select(orthogroup_id, species, participation, og_type) %>%
  pivot_wider(names_from = species, values_from = participation) %>%
  drop_na(Scil, Slac)
mtm_wide_participation

mtm_wide_redundancy <- og_breadth_metrics %>%
  filter(og_type != "1 to 1") %>%
  select(orthogroup_id, species, redundancy, og_type) %>%
  pivot_wider(names_from = species, values_from = redundancy) %>%
  drop_na(Scil, Slac)
mtm_wide_redundancy

paired_tests_mtm <- tibble(
  metric = c(rep("participation", 3), rep("redundancy", 3)),
  og_type = rep(c("Single copy\nin Sycon", "Single copy\nin Spongilla", "Many to many"), 2),
  # n_pairs = c(nrow(mtm_wide_participation),
  #             nrow(mtm_wide_redundancy)),
  p_value = c(wilcox.test(mtm_wide_participation %>% filter(og_type == "Single copy\nin Sycon") %>%
                            pull(Scil),
                          mtm_wide_participation %>% filter(og_type == "Single copy\nin Sycon") %>%
                            pull(Slac), paired = TRUE)$p.value,
              wilcox.test(mtm_wide_participation %>% filter(og_type == "Single copy\nin Spongilla") %>%
                            pull(Scil),
                          mtm_wide_participation %>% filter(og_type == "Single copy\nin Spongilla") %>%
                            pull(Slac), paired = TRUE)$p.value,
              wilcox.test(mtm_wide_participation %>% filter(og_type == "Many to many") %>%
                            pull(Scil),
                          mtm_wide_participation %>% filter(og_type == "Many to many") %>%
                            pull(Slac), paired = TRUE)$p.value,
              wilcox.test(mtm_wide_redundancy %>% filter(og_type == "Single copy\nin Sycon") %>%
                            pull(Scil),
                          mtm_wide_redundancy %>% filter(og_type == "Single copy\nin Sycon") %>%
                            pull(Slac),
                          paired = TRUE)$p.value,
              wilcox.test(mtm_wide_redundancy %>% filter(og_type == "Single copy\nin Spongilla") %>%
                            pull(Scil),
                          mtm_wide_redundancy %>% filter(og_type == "Single copy\nin Spongilla") %>%
                            pull(Slac),
                          paired = TRUE)$p.value,
              wilcox.test(mtm_wide_redundancy %>% filter(og_type == "Many to many") %>%
                            pull(Scil),
                          mtm_wide_redundancy %>% filter(og_type == "Many to many") %>%
                            pull(Slac),
                          paired = TRUE)$p.value)) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"),
         label = case_when(p_adj < 0.001 ~ "***",
                           p_adj < 0.01 ~ "**",
                           p_adj < 0.05 ~ "*",
                           TRUE ~ "ns"))
paired_tests_mtm

mtm_long <- bind_rows(
  mtm_wide_participation %>% pivot_longer(c(Scil, Slac), names_to = "species", values_to = "value") %>% mutate(metric = "participation"),
  mtm_wide_redundancy    %>% pivot_longer(c(Scil, Slac), names_to = "species", values_to = "value") %>% mutate(metric = "redundancy")) %>%
  mutate(species_label = if_else(species == "Scil", "S. ciliatum", "S. lacustris"),
         og_type = factor(og_type, levels = c("Single copy\nin Sycon",
                                              "Single copy\nin Spongilla",
                                              "Many to many")))
mtm_plot <- mtm_long %>%
  ggplot(aes(species_label, value)) +
  geom_line(aes(group = orthogroup_id), alpha = 0.4, colour = "grey70") +
  geom_jitter(aes(col = species_label), width = 0.15, height = 0.01) +
  geom_boxplot(outliers = FALSE, fill = alpha("white", 0),
               linewidth = 0.7, width = 0.4) +
  
  scale_color_manual(values = c("S. ciliatum" = "#ff9ebb", "S. lacustris" = "#b9375e"), guide = "none") +
  
  facet_grid(metric ~ og_type, scales = "free_y",
             labeller = labeller(og_type = as_labeller(c("Single copy\nin Sycon" = paste0("Single copy in Sycon\n(n = ", count_og_type[["Single copy\nin Sycon"]], ")"),
                                                         "Single copy\nin Spongilla" = paste0("Single copy in Spongilla\n(n = ", count_og_type[["Single copy\nin Spongilla"]], ")"),
                                                         "Many to many" = paste0("Many to many\n(n = ", count_og_type[["Many to many"]], ")"))),
                                 metric = as_labeller(c("participation" = "Participation",
                                                        "redundancy" = "Redundancy")))) +
  geom_text(data = paired_tests_mtm,
            aes(x = 1.5, y = c(rep(1.2, 3), rep(2.2, 3)),
                label = label), vjust = 1,
            inherit.aes = FALSE, size = 4) +
  
  theme_bw(base_size = 12) +
  theme_for_plots +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, face = "italic", hjust = 1))
mtm_plot

ggsave("17_comparisons_1to1_orthologues/mtm_boxplot_statistics.png",
       mtm_plot, device = "png",
       height = 6.5, width = 6, dpi = 300, unit = "in", bg = "white")
ggsave("17_comparisons_1to1_orthologues/mtm_boxplot_statistics.pdf",
       mtm_plot, device = cairo_pdf,
       height = 6.5, width = 6, dpi = 300, unit = "in", bg = "white")



#############################################################
#     BREADTH VS. N_GENES REGRESSION                        #
#############################################################

orthogroup_full_set_with_geneCounts_clusterCounts %>%
  distinct(orthogroup_id, species, species_label, n_genes, n_cluster_prop) %>%
  ggplot(aes(n_genes, n_cluster_prop, fill = species_label)) +
  geom_jitter(aes(col = species_label), width = 0.1, height = 0.1) +
  scale_color_manual(values = c("S. ciliatum" = alpha("#ff9ebb", 0.7), "S. lacustris" = alpha("#b9375e", 0.7))) +
  scale_fill_manual(values = c("S. ciliatum" = alpha("#ff9ebb", 0.7), "S. lacustris" = alpha("#b9375e", 0.7))) +
  ggnewscale::new_scale_color() +
  geom_smooth(aes(col = species_label), method = "lm") +
  scale_color_manual(values = c("S. ciliatum" = "#ff9ebb", "S. lacustris" = "#b9375e"), guide = "none") +
  labs(y = "Proportion of clusters", x = "Number of genes") +
  theme_bw() +
  theme_for_plots


##############################################
#     CLUSTER TRANSCRIPTIONAL UNIQUENESS     #
##############################################

# compute the number of exclusive 1 to 1 OGs for each cluster
cluster_specific_1to1 <- orthogroup_full_set_with_clusterCount %>%
  group_by(orthogroup_id) %>%
  # remove all OGs where at least one gene is expressed in any cluster
  filter(!any(n_clusters == 0)) %>%
  ungroup() %>%
  
  left_join(og_type_table, by = "orthogroup_id") %>%
  
  # calculate the total number of markers per cluster
  group_by(cluster) %>%
  mutate(total_markers_perCluster = n_distinct(orthogroup_id)) %>%
  ungroup() %>%
  
  filter(og_type == "1 to 1") %>%
  
  # keep only orthogroup/species combos that mark exactly 1 cluster
  group_by(orthogroup_id, species) %>%
  filter(n_distinct(cluster) == 1) %>%
  ungroup() %>%
  
  group_by(cluster) %>%
  mutate(n_1to1_specific_perCluster = n_distinct(orthogroup_id),
         prop_1to1_specific_perCluster = n_1to1_specific_perCluster/total_markers_perCluster) %>%
  ungroup() %>%
  distinct(cluster, prop_1to1_specific_perCluster, n_1to1_specific_perCluster, total_markers_perCluster)
cluster_specific_1to1

# subset the integrated object to only scil
scil_samap <- ScilSlac_samap %>%
  subset(species == "Scil")

# update metadata with the proportion of exclusive 1-to-1 OGs
scil_samap[[]] <- scil_samap[[]] %>%
  rownames_to_column() %>%
  mutate(Scil_leiden_clusters = paste0("Scil_", Scil_leiden_clusters)) %>%
  left_join(cluster_specific_1to1,
            by = join_by("Scil_leiden_clusters" == "cluster")) %>%
  replace(is.na(.), 0) %>%
  column_to_rownames()
scil_samap[[]] %>% View


# subset the integrated object to only slac
slac_samap <- ScilSlac_samap %>%
  subset(species == "Slac")

# update metadata with the proportion of exclusive 1-to-1 OGs
slac_samap[[]] <- slac_samap[[]] %>%
  rownames_to_column() %>%
  mutate(Slac_leiden_clusters = paste0("Slac_", Slac_leiden_clusters)) %>%
  left_join(cluster_specific_1to1,
            by = join_by("Slac_leiden_clusters" == "cluster")) %>%
  replace(is.na(.), 0) %>%
  column_to_rownames()
slac_samap[[]]

# plot UMAPs
panel_umap <- ScilSlac_samap %>%
  DimPlot2(group.by = "species", pt.size = 0.8, cols = c("#ff9ebb", "#b9375e"),
           theme = theme_umap_arrows()) +
  scil_samap %>%
  DimPlot2(features = "prop_1to1_specific_perCluster", pt.size = 0.8,
           theme = theme_umap_arrows()) +
  slac_samap %>%
  DimPlot2(features = "prop_1to1_specific_perCluster", pt.size = 0.8,
           theme = theme_umap_arrows())

panel_umap

ggsave("17_comparisons_1to1_orthologues/panel_cluster_uniqueness.png",
       panel_umap, device = "png",
       height = 3*1.5, width = 10*1.5, dpi = 300, unit = "in", bg = "white")
ggsave("17_comparisons_1to1_orthologues/panel_cluster_uniqueness.pdf",
       panel_umap, device = cairo_pdf,
       height = 3*1.5, width = 10*1.5, dpi = 300, unit = "in", bg = "white")

scilslac_umap_integrated_raw <- DimPlot(ScilSlac_samap, group.by = "species")
scil_umap_prop_1to1_raw <- FeaturePlot(scil_samap, features = "prop_1to1_specific_perCluster")
slac_umap_prop_1to1_raw <- FeaturePlot(slac_samap, features = "prop_1to1_specific_perCluster")

scilslac_umap_integrated <- scilslac_umap_integrated_raw@data %>%
  arrange(desc(species)) %>%
  mutate(species = case_when(species == "Scil" ~ "S. ciliatum",
                             species == "Slac" ~ "S. lacustris",
                             TRUE ~ species)) %>%
  
  ggplot(aes(Xumap_1, Xumap_2, col = species)) +
  geom_point(size = 0.01) +
  
  scale_color_manual(values = c("#ff9ebb", "#b9375e")) +
  
  guides(color = guide_legend(label.position = "bottom",
                              override.aes = list(size = 4))) +
  
  labs(title = "SAMap integrated space", col = "Species") +
  
  umap_arrows +
  theme_bw(base_size = 12) +
  theme_for_UMAPS +
  theme(legend.position = "bottom",
        legend.title = element_text(hjust = 0.5, size = 10),
        legend.title.position = "top",
        legend.text = element_text(face = "italic"))
scilslac_umap_integrated

# scil_umap_prop_1to1 <- scil_umap_prop_1to1_raw@data %>%
#   arrange(prop_1to1_specific_perCluster) %>%
#   
#   ggplot(aes(Xumap_1, Xumap_2, col = prop_1to1_specific_perCluster)) +
#   geom_point(size = 0.8) +
#   
#   labs(title = expression(paste(bolditalic("S. ciliatum"), bold(" cell subset"))),
#        col = "Transcriptional uniqueness\n(1-to-1 orthologues)") +
#   
#   scale_color_distiller(palette = "Blues", direction = 1, breaks = c(0, 0.25, 0.5)) +
#   guides(size = guide_legend(label.position = "bottom"),
#          color = guide_colorbar(barheight = 1)) +
#   
#   theme_bw(base_size = 12) +
#   theme_for_UMAPS +
#   theme(legend.position = "bottom",
#         legend.title = element_text(hjust = 0.5, size = 10),
#         legend.title.position = "top")
# scil_umap_prop_1to1
# 
# 
# slac_umap_prop_1to1 <- slac_umap_prop_1to1_raw@data %>%
#   arrange(prop_1to1_specific_perCluster) %>%
#   
#   ggplot(aes(Xumap_1, Xumap_2, col = prop_1to1_specific_perCluster)) +
#   geom_point(size = 0.8) +
#   
#   labs(title = expression(paste(bolditalic("S. lacustris"), bold(" cell subset"))),
#        col = "Transcriptional uniqueness\n(1-to-1 orthologues)") +
#   
#   scale_color_distiller(palette = "Blues", direction = 1, breaks = c(0, 0.03, 0.06)) +
#   guides(size = guide_legend(label.position = "bottom"),
#          color = guide_colorbar(barheight = 1)) +
#   
#   theme_bw(base_size = 12) +
#   theme_for_UMAPS +
#   theme(legend.position = "bottom",
#         legend.title = element_text(hjust = 0.5, size = 10),
#         legend.title.position = "top")
# slac_umap_prop_1to1

scilslac_umap_prop_1to1_faceted <- bind_rows(slac_umap_prop_1to1_raw@data %>%
                                               mutate(species = "Slac"),
                                             scil_umap_prop_1to1_raw@data %>%
                                               mutate(species = "Scil")) %>%
  arrange(prop_1to1_specific_perCluster) %>%
  
  ggplot(aes(Xumap_1, Xumap_2, col = prop_1to1_specific_perCluster)) +
  geom_point(size = 0.01) +
  
  labs(title = expression(paste(bolditalic("S. lacustris"), bold(" cell subset"))),
       col = "TUi") +
  
  scale_color_distiller(palette = "Blues", direction = 1,
                        breaks = c(0, 0.125, 0.25)) +
  
  guides(size = guide_legend(label.position = "bottom"),
         # color = guide_colorbar(barheight = 1)
         ) +
  
  facet_wrap(~ species,
             labeller = labeller(species = as_labeller(
               c("Scil" = "S. ciliatum cells",
                 "Slac" = "S. lacustris cells")))) +
  
  theme_bw(base_size = 12) +
  theme_for_UMAPS +
  theme(strip.text = element_text(color = "black", face = "bold", hjust = 0.5),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.clip = "off",
        plot.title = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(hjust = 0.5, size = 10),
        legend.title.position = "top")
scilslac_umap_prop_1to1_faceted

panel_umap <- ggpubr::ggarrange(scilslac_umap_integrated, scilslac_umap_prop_1to1_faceted,
                                widths = c(1,2))
panel_umap

ggsave("17_comparisons_1to1_orthologues/panel_cluster_uniqueness.png",
       panel_umap, device = "png",
       height = 4, width = 8, dpi = 300, unit = "in", bg = "white")
ggsave("17_comparisons_1to1_orthologues/panel_cluster_uniqueness.pdf",
       panel_umap, device = cairo_pdf,
       height = 4, width = 8, dpi = 300, unit = "in", bg = "white")


###########################################
#     CORRELATION WITH MAPPING SCORES     #
###########################################

mappScore_per_tui <- read.table("05NEW_SAMap_porifera/01_mapping_scores/ScilSlac_leiden3Clusters_leiden_100topCells_samapMappingTable.tsv",
           header = TRUE, sep = "\t") %>%
  
  # create a long-format dataframe
  pivot_longer(-X, names_to = "target_ID", values_to = "mapp_score") %>%
  rename("source_ID" = "X") %>% 
  
  # filter out rows with mapping scores below the threshold
  filter(mapp_score > 0) %>%
  
  group_by(source_ID) %>%
  slice_max(order_by = mapp_score) %>%
  ungroup() %>%
  
  left_join(cluster_specific_1to1,
            by = join_by("source_ID" == "cluster")) %>%
  replace_na(list(prop_1to1_specific_perCluster = 0,
                  n_1to1_specific_perCluster = 0,
                  total_markers_perCluster = 0)) %>%
  separate(source_ID, into = c("source_species"), sep = "_",
           remove = FALSE)

lm_plots <- mappScore_per_tui %>%
  ggplot(aes(mapp_score, prop_1to1_specific_perCluster,
             fill = source_species, col = source_species)) +
  
  geom_point() +
  geom_smooth(method = "lm", alpha = 0.2) +
  
  labs(x = "SAMap mapping score", y = "TUi") +
  scale_color_manual(values = c("#ff9ebb", "#b9375e")) +
  scale_fill_manual(values = c("#ff9ebb", "#b9375e")) +
  
  scale_x_continuous(limits = c(0, 1),
                     breaks = c(0, 0.5, 1)) +
  
  scale_y_continuous(breaks = function(x) {
    c(0, max(x, na.rm = TRUE) / 2, max(x, na.rm = TRUE))},
    labels = scales::number_format(accuracy = 0.01)) +

  facet_wrap(~ source_species, scale = "free_y", nrow = 2) +

  theme_bw(base_size = 12) +
  theme(plot.background = element_rect(fill = "transparent", colour = NA),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey90", lineend = "round"),
        panel.border = element_rect(colour = "black", linewidth = .6),
        legend.position = "none",
        axis.line = element_blank(),
        axis.ticks = element_line(colour = "black", linewidth = .4),
        axis.ticks.length = unit(0.10, "cm"),
        axis.text.x = element_text(color = "black", margin = margin(t = 4, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(color = "black", margin = margin(t = 0, r = 4, b = 0, l = 0)),
        axis.title.y = element_text(angle = 90, size = 12, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(angle = 0, size = 12, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        strip.text = element_blank(),
        strip.background = element_blank())
lm_plots

lm <- mappScore_per_tui %>%
  group_by(source_species) %>%
  do(model = lm(prop_1to1_specific_perCluster ~ mapp_score, data = .))

ggsave("17_comparisons_1to1_orthologues/linear_models_plots.png",
       lm_plots, device = "png",
       height = 4, width = 3, dpi = 300, unit = "in", bg = "white")
ggsave("17_comparisons_1to1_orthologues/linear_models_plots.pdf",
       lm_plots, device = cairo_pdf,
       height = 1, width = 2, dpi = 300, unit = "in", bg = "white")


############################
#     CHECK CLUSTER ID     #
############################

scil_samap[[]] %>%
  select(Scil_leiden_clusters, Scil_seurat_clusters) %>%
  left_join(read.table("07_notableGenes_clusterAnnotation/cluster_identity.tsv", header = TRUE) %>%
              mutate(cluster_ID = as.factor(cluster_ID)),
            by = join_by("Scil_seurat_clusters" == "cluster_ID")) %>%
  count(Scil_leiden_clusters, Scil_seurat_clusters, value) %>%
  arrange(Scil_leiden_clusters, desc(n)) %>% 
  filter(Scil_leiden_clusters %in% paste0("Scil_", c(24,34,26,8,178,28,33))) %>%
  View

DimPlot2(scil_samap, features = "Scil_leiden_clusters",
         label = TRUE, repel = TRUE, theme = NoLegend())
DimPlot2(slac_samap, features = "Slac_leiden_clusters",
         label = TRUE, repel = TRUE, theme = NoLegend())

slac_samap[[]] %>%
  select(Slac_leiden_clusters, Slac_seurat_clusters_1, Slac_seurat_clusters_2) %>%
  left_join(
    data.frame(
      cell_type = c("Slac_0$" = "Archaeocytes S1", "Slac_1$" = "Archaeocytes S2",
                    "Slac_2$" = "Archaeocytes S3", "Slac_3$" = "Archaeocytes S4",
                    "Slac_4$" = "Archaeocytes S5", "Slac_5$" = "Archaeocytes S4",
                    "Slac_9$" = "Myopeptidocytes", "Slac_14$" = "Metabolocytes",
                    "Slac_15$" = "Myopeptidocytes", "Slac_17$" = "Pinacocytes",
                    "Slac_21$" = "Archaeocytes S6", "Slac_23$" = "Archaeocytes-like",
                    "Slac_24$" = "Myopeptidocytes", "Slac_25$" = "Choanocytes/-blasts",
                    "Slac_30$" = "Sclerocytes", "Slac_31$" = "Amoebocytes/Neuroid",
                    "Slac_35$" = "Basopinacocytes", "Slac_36$" = "Mesocytes",
                    "Slac_37$" = "Granulocytes-like", "Slac_41$" = "Mesocytes")
    ) %>%
      rownames_to_column() %>%
      mutate(rowname = str_remove(rowname, "\\$"),
             rowname = str_remove(rowname, "Slac_")),
    by = join_by(Slac_seurat_clusters_2 == rowname)
  ) %>%
  select(-Slac_seurat_clusters_1) %>% 
  count(Slac_leiden_clusters, Slac_seurat_clusters_2, cell_type) %>%
  arrange(Slac_leiden_clusters, desc(n)) %>%
  filter(Slac_leiden_clusters %in% paste0("Slac_", c(18, 19, 13, 15,21, 6))) %>%
  View
