#!/usr/bin/env Rscript

setwd("/data/evassvis/fn76/sycon/sycon_clusterAnnotation/sycon_cluster_analyses/")

library(tidyverse)
library(Seurat)
library(SeuratExtend)
library(VennDiagram)


###############################
#     DEFINE PLOT THEMES      #
###############################

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
  strip.placement = "outside",
  strip.background = element_blank(),
  strip.clip = "off"
)

umap_arrows <- list(
  
  annotation_custom(grob = grid::segmentsGrob(x0 = unit(0, "mm"), x1 = unit(12, "mm"),
                                              y0 = unit(0, "mm"), y1 = unit(0, "mm"),
                                              arrow = arrow(length = unit(2.5, "mm"),
                                                            ends = "last", type = "open"),
                                              gp = grid::gpar(col = "black", fill = "black", lwd = 1))),
    
    annotation_custom(grob = grid::segmentsGrob(x0 = unit(0, "mm"), x1 = unit(0, "mm"),
                                                y0 = unit(0, "mm"), y1 = unit(12, "mm"),
                                                arrow = arrow(length = unit(2.5, "mm"),
                                                              ends = "last", type = "open"),
                                                gp = grid::gpar(col = "black", fill = "black", lwd = 1))),
    
    annotation_custom(grob = grid::textGrob(label = "UMAP 1",
                                            x = unit(0, "mm"), y = unit(0, "mm") - unit(2.5, "mm"),
                                            just = c(0, 1), gp = grid::gpar(fontsize = 10))),
    
    annotation_custom(grob = grid::textGrob(label = "UMAP 2",
                                            x = unit(0, "mm") - unit(2.5, "mm"), y = unit(0, "mm"),
                                            just = c(0, 0), rot = 90, gp = grid::gpar(fontsize = 10))),
    
    coord_cartesian(clip = "off"))


theme_for_UMAPS <- theme(
  # aspect.ratio = 1,
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

# scilslac <- schard::h5ad2seurat("05NEW_SAMap_porifera/ScilSlac_leiden3Clusters_samap.h5ad")
# aquescilslac <- schard::h5ad2seurat("05NEW_SAMap_porifera/AqueScilSlac_leiden3Clusters_samap.h5ad")
# aqueslac <- schard::h5ad2seurat("05NEW_SAMap_porifera/AqueSlac_leiden3Clusters_samap.h5ad")

# import orthogroup data
PCA_1to1 <- read.table("17_comparisons_1to1_orthologues/01_orthologue_tables/1_to_1_Sycon_Spongilla.tsv",
                       header = TRUE, sep = "\t", na.strings = c("", ".")) %>%
  select(-c(required_nodes_present, required_single_nodes_present, collected_specs_present)) %>%
  mutate(level = "PCA",
         # og_type = "1-to-1"
         ) %>%
  as_tibble()
PCA_1to1

MCA_1to1 <- read.table("17_comparisons_1to1_orthologues/01_orthologue_tables/1_to_1_Sycon_Spongilla_present_on_the_metazoan_stem.tsv",
                       header = TRUE, sep = "\t", na.strings = c("", ".")) %>%
  select(-c(required_nodes_present, required_single_nodes_present, collected_specs_present)) %>%
  mutate(level = "MCA",
         # og_type = "1-to-1"
         ) %>%
  as_tibble()
MCA_1to1


PCA_1toMany <- read.table("17_comparisons_1to1_orthologues/01_orthologue_tables/1_to_many_Sycon_Spongilla.tsv",
                          header = TRUE, sep = "\t", na.strings = c("", ".")) %>%
  select(-c(required_nodes_present, required_single_nodes_present, collected_specs_present)) %>%
  mutate(level = "PCA",
         # og_type = "1-to-many"
         ) %>%
  as_tibble()
PCA_1toMany

MCA_1toMany <- read.table("17_comparisons_1to1_orthologues/01_orthologue_tables/1_to_many_Sycon_Spongilla_present_on_the_metazoan_stem.tsv",
                          header = TRUE, sep = "\t", na.strings = c("", ".")) %>%
  select(-c(required_nodes_present, required_single_nodes_present, collected_specs_present)) %>%
  mutate(level = "MCA",
         # og_type = "1-to-many"
         ) %>%
  as_tibble()
MCA_1toMany

# check how they intersect
venn.diagram(list("PCA_1to1" = PCA_1to1$orthogroup_id,
                  "MCA_1to1" = MCA_1to1$orthogroup_id,
                  "PCA_1toMany" = PCA_1toMany$orthogroup_id,
                  "MCA_1toMany" = MCA_1toMany$orthogroup_id),
             filename = NULL, fill = c("red", "green", "blue", "yellow"))

# load slac gene name conversion table
# this is because gene names and protein names do not correspond of course
# so e.g. the transcript and protein id of the gene ENSLPGG00000000002 are ENSLPGT00000000015 and ENSLPGP00000000013, respectively
# I wonder what the hell people are thinking when they write programs...
slac_gene_names_conversion <- read.table("../spongilla_remapping/00_input/slac_genome/genes_to_transcript.tsv",
                                         sep = "\t", na.strings = "N/A") %>%
  left_join(read.table("../spongilla_remapping/00_input/slac_genome/transcripts_to_cdss.tsv",
                       sep = "\t", na.strings = "N/A"),
            by = join_by("V1" == "V2")) %>%
  as_tibble() %>%
  rename("transcript_id" = "V1",
         "gene_id" = "V2",
         "protein_id" = "V1.y",) %>%
  arrange(gene_id)
slac_gene_names_conversion

# # load SAMap connecting gene pairs table
# ScilSlac_samapGenePairs <- read.table("05_SAMap_porifera/02_gene_pairs/02_threshold04_leidenClusters/ScilSlac_leiden3Clusters_all_samapGenePairs.tsv",
#            fill = TRUE, na.strings = "", header = TRUE) %>%
#   rename_with(~ paste0(., "_gene_pair"), .cols = !matches("_pval[12]$")) %>%
#   pivot_longer(cols = everything(),
#                names_to = c("mapping_clusters", ".value"),
#                names_pattern = "^(.*)_(gene_pair|pval1|pval2)$") %>%
#   drop_na()


################################################
#     COUNT HOW MANY GENES PER ORTHOGROUPS     #
################################################

# get the full set of orthogroups, each gene on a different row
orthogroup_full_set <- bind_rows(MCA_1to1, MCA_1toMany, PCA_1to1, PCA_1toMany) %>%
  distinct(orthogroup_id, .keep_all = TRUE) %>%
  pivot_longer(c(LDemCA__Spongilla_lacustris, LCalcCA__Sycon_ciliatum),
               names_to = "species", values_to = "gene_ID") %>%
  mutate(species = case_when(str_detect(species, "Spongilla") ~ "Slac",
                             str_detect(species, "Sycon") ~ "Scil",
                             TRUE ~ species)) %>%
  # separate each paralogue/isoform on a different row
  separate_rows(gene_ID, sep = ";") %>%
  
  # remove spongilla ".1$"
  mutate(gene_ID = str_remove(gene_ID, "\\.[0-9]+$")) %>%
  
  # add information on gene ids for spongilla, which are used in the single cell data
  left_join(slac_gene_names_conversion,
            by = join_by("gene_ID" == "protein_id")) %>%
  mutate(gene_id = case_when(is.na(gene_id) ~ gene_ID,
                             TRUE ~ gene_id)) %>%
  select(-transcript_id) %>%
  rename("protein_id" = "gene_ID")
orthogroup_full_set

# the code from here was modified because in Spongilla, paralogues are sometimes actually isoforms
# so we need to keep both protein_id and gene_id to be able to compure orthogroup type

# the main thing that is changing is the gene count per orthogroup,
# cos of course isoforms do not contribute to gene count
# so we need to convert from protein ids to gene ids and count those
orthogroup_full_set_with_geneCounts <- orthogroup_full_set %>%
  group_by(orthogroup_id, species) %>%
  mutate(n_genes = n_distinct(gene_id))

orthogroup_full_set_with_geneCounts


#######################################
#     LOAD AND PROCESS SAMap DATA     #
#######################################

# SKIP THIS SNIPPET IF RE-RUNNING

# load the integrated SAMap object
ScilSlac_samap <- sceasy::convertFormat("05_SAMap_porifera/ScilSlac_leiden3Clusters_samap.h5ad",
                                        from = "anndata", to = "seurat")

tmp <- ScilSlac_samap %>% DimPlot(group.by = "Scil_orig.ident")

tmp@data %>%
  filter(Scil_orig.ident != "unassigned") %>%
  ggplot(aes(UMAP_1, UMAP_2, col = Scil_orig.ident)) +
  geom_point(size = 0.8)

# plot and save the integrated UMAP
ScilSlac_samap %>%
  SeuratExtend::DimPlot2(group.by = "species", cols = "bright",
                         theme = theme_umap_arrows()) %>%
  ggsave(file = "17_comparisons_1to1_orthologues/SAMap_integrated_umap.png",
         device = "png", width = 8, height = 6, dpi = 300, unit = "in", bg = "white")

ScilSlac_samap %>%
  SeuratExtend::DimPlot2(group.by = "species", cols = "bright",
                         theme = theme_umap_arrows()) %>%
  ggsave(file = "17_comparisons_1to1_orthologues/SAMap_integrated_umap.pdf",
         device = cairo_pdf, width = 8, height = 6, dpi = 300, unit = "in", bg = "white")

# plot only scil with leiden clusters in the integrated space
ScilSlac_samap %>%
  subset(species == "Scil") %>%
  SeuratExtend::DimPlot2(group.by = "Scil_leiden_clusters", cols = "bright",
                         theme = theme_umap_arrows()) %>%
  ggsave(file = "17_comparisons_1to1_orthologues/SAMap_scil_leidenClusters.png",
         device = "png", width = 9, height = 6, dpi = 300, unit = "in", bg = "white")

ScilSlac_samap %>%
  subset(species == "Scil") %>%
  SeuratExtend::DimPlot2(group.by = "Scil_leiden_clusters", cols = "bright",
                         theme = theme_umap_arrows()) %>%
  ggsave(file = "17_comparisons_1to1_orthologues/SAMap_scil_leidenClusters.pdf",
         device = cairo_pdf, width = 9, height = 6, dpi = 300, unit = "in", bg = "white")

# plot only slac with leiden clusters in the integrated space
ScilSlac_samap %>%
  subset(species == "Slac") %>%
  SeuratExtend::DimPlot2(group.by = "Slac_leiden_clusters", cols = "bright",
                         theme = theme_umap_arrows()) %>%
  ggsave(file = "17_comparisons_1to1_orthologues/SAMap_slac_leidenClusters.png",
         device = "png", width = 9, height = 6, dpi = 300, unit = "in", bg = "white")

ScilSlac_samap %>%
  subset(species == "Slac") %>%
  SeuratExtend::DimPlot2(group.by = "Slac_leiden_clusters", cols = "bright",
                         theme = theme_umap_arrows()) %>%
  ggsave(file = "17_comparisons_1to1_orthologues/SAMap_slac_leidenClusters.pdf",
         device = cairo_pdf, width = 9, height = 6, dpi = 300, unit = "in", bg = "white")


# subset the integrated object to only scil
scil_samap <- ScilSlac_samap %>%
  subset(species == "Scil") %>%
  SCTransform()

# find markers of scil leiden clusters
scil_samap_allMarkers <- scil_samap %>%
  FindAllMarkers(group.by = "Scil_leiden_clusters") %>%
  rename("leiden_clusters" = "cluster") %>%
  mutate(gene = str_replace(gene, "Scil-", ""))

write.table(scil_samap_allMarkers, "17_comparisons_1to1_orthologues/scil_samap_allMarkers_leidenClusters.tsv",
            col.names = TRUE, sep = "\t", quote = FALSE, row.names = FALSE)

# subset the integrated object to only slac
slac_samap <- ScilSlac_samap %>%
  subset(species == "Slac") %>%
  SCTransform()

# find markers of slac leiden clusters
slac_samap_allMarkers <- slac_samap %>%
  FindAllMarkers(group.by = "Slac_leiden_clusters") %>%
  rename("leiden_clusters" = "cluster") %>%
  mutate(gene = str_replace(gene, "Slac-", "")) %>%
  left_join(slac_gene_names_conversion, by = join_by("gene" == "gene_id"),
            relationship = "many-to-many")

write.table(slac_samap_allMarkers, "17_comparisons_1to1_orthologues/slac_samap_allMarkers_leidenClusters.tsv",
            col.names = TRUE, sep = "\t", quote = FALSE, row.names = FALSE)


#####################################################################
#     COMPUTE THE NUMBER OF CLUSTERS WHERE EACH OG IS EXPRESSED     #
#####################################################################

# load all markers tables
scil_samap_allMarkers <- read.table("17_comparisons_1to1_orthologues/scil_samap_allMarkers_leidenClusters.tsv",
                                    header = TRUE, sep = "\t") %>% tibble()
slac_samap_allMarkers <- read.table("17_comparisons_1to1_orthologues/slac_samap_allMarkers_leidenClusters.tsv",
                                    header = TRUE, sep = "\t") %>% tibble()

# compile a unique table for both slac and spongilla
allMarkers_per_cluster <- scil_samap_allMarkers %>%
  
  # add the species ID before the cluster name
  mutate(leiden_clusters = paste0("Scil_", leiden_clusters),
         protein_id = gene) %>%
  
  # append slac markers
  add_row(slac_samap_allMarkers %>%
            # keep only the gene_id column
            select(-transcript_id) %>%
            mutate(leiden_clusters = paste0("Slac_", leiden_clusters)), ) %>%
  # filter out non significant markers and those with log2FC < 1
  filter(p_val_adj < 0.05,
         avg_log2FC > 1) %>%
  select(leiden_clusters, gene, protein_id) %>%
  rename("gene_id" = "gene")
allMarkers_per_cluster


# considering the above about Slac isoforms, the cluster count per orthogroup stays the same,
# cos regardless of paralogues being true or isoforms, still they are expressed in clusters
# so in theory nothing changes from the cluster count perspective

# add indication of cluster expression
orthogroup_full_set_withClusterExpression <- orthogroup_full_set %>%
  
  # add cluster expression based on compiled markers
  left_join(allMarkers_per_cluster,
            by = join_by(gene_id, protein_id),
            relationship = "many-to-many") %>%
  group_by(orthogroup_id, species) %>%
  
  # count in how many clusters each og in each species is expressed
  mutate(n_clusters = n_distinct(leiden_clusters, na.rm = TRUE)) %>%
  ungroup()
orthogroup_full_set_withClusterExpression

# tidy up data frame
orthogroup_full_set_with_geneCounts_clusterCounts <- orthogroup_full_set_with_geneCounts %>%
  
  # add cluster numbers to the tibble with gene numbers
  left_join(orthogroup_full_set_withClusterExpression %>%
              # remove reduntant rows on cluster expression
              select(orthogroup_id, species, n_clusters) %>%
              distinct(),
            by = join_by(orthogroup_id, species)) %>%
  
  # add information on orthogroup type
  left_join((.) %>%
              select(orthogroup_id, n_genes, species) %>%
              distinct() %>%
              pivot_wider(names_from = species, values_from = n_genes) %>%
              mutate(og_type = case_when(Scil == 1 & Slac == 1 ~ "1 to 1",
                                         Scil > 1 & Slac == 1 ~ "Single copy\nin Spongilla",
                                         Scil == 1 & Slac > 1 ~ "Single copy\nin Sycon",
                                         TRUE ~ "Many to many")) %>%
              select(orthogroup_id, og_type),
            by = join_by(orthogroup_id)) %>%
  
  # compute the proportion of clusters 
  mutate(n_cluster_prop = case_when(species == "Slac" ~ n_clusters/33,
                                    species == "Scil" ~ n_clusters/42)) %>%
  
  # remove orthogroups with genes that are not expressed in any cluster
  group_by(orthogroup_id) %>%
  filter(any(n_clusters > 0)) %>%
  ungroup() %>%
  
  # edit species name and order og_type factor
  mutate(species = case_when(species == "Scil" ~ "S. ciliatum",
                             species == "Slac" ~ "S. lacustris",
                             TRUE ~ species),
         og_type = factor(og_type, levels = c("1 to 1", "Single copy\nin Sycon",
                                              "Single copy\nin Spongilla", "Many to many")))
  # filter(orthogroup_id %in% c("OG0002796", "OG0009672"))
  
orthogroup_full_set_with_geneCounts_clusterCounts


#############################################################
#     BOXPLOTS WITH PCA ORTHOGROUPS WITH EITHER SPECIES     #
#############################################################

# compute orthogroups with proportion of clusters
boxplots_per_og_type <- orthogroup_full_set_with_geneCounts_clusterCounts %>%
  
  ggplot(aes(x = species, y = n_cluster_prop)) +
  
  geom_jitter(aes(col = species), width = 0.2, height = 0.01) +
  geom_boxplot(outliers = FALSE, fill = alpha("white", 0.8),
               width = 0.8) +
  
  scale_color_manual(values = c("#ff9ebb", "#b9375e"), guide = "none") +
  
  labs(y = "Proportion of clusters") +
  
  ggpubr::stat_compare_means(aes(group = species), label = "p.signif") +
  
  facet_wrap(~og_type, ncol = 4) +
  theme_bw(base_size = 12) +
  theme_for_plots +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, face = "italic", hjust = 1))
boxplots_per_og_type

# compute distribution of proportion of clusters
ridges_per_species <- orthogroup_full_set_with_geneCounts_clusterCounts %>%
  
  ggplot(aes(x = n_cluster_prop, y = species, fill = species, col = species)) +
  ggridges::geom_density_ridges(bandwidth = 0.02,
                                alpha = 0.8, linewidth = 0.6) +
  
  scale_color_manual(values = c("#ff9ebb", "#b9375e"), guide = "none") +
  scale_fill_manual(values = c("#ff9ebb", "#b9375e"), guide = "none") +
  scale_x_continuous(breaks = seq(0, 0.3, 0.1),
                     limits = c(0, 0.35)) +
  
  coord_flip() +
  
  facet_wrap(~ species) +

  theme_bw(base_size = 12) +
  theme_for_plots +
  theme(panel.border = element_blank(),
        # panel.grid.major.y = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        # axis.text.y = element_blank(),
        # axis.text.x = element_text(angle = 45, face = "italic", hjust = 1),
        axis.text.x = element_blank(),
        # axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank())
ridges_per_species

panel <- ( boxplots_per_og_type + ridges_per_species ) +
  plot_layout(axes = "collect",
              widths = c(1, 0.5))
panel

ggsave("17_comparisons_1to1_orthologues/boxplot_orthogroup_categories.png",
       panel, device = "png",
       height = 4*1.5, width = 10*1.5, dpi = 300, unit = "in", bg = "white")
ggsave("17_comparisons_1to1_orthologues/boxplot_orthogroup_categories.pdf",
       panel, device = cairo_pdf,
       height = 4*1.5, width = 10*1.5, dpi = 300, unit = "in", bg = "white")

# compute linear regressions
orthogroup_full_set_with_geneCounts_clusterCounts %>%
  
  ggplot(aes(n_genes, n_cluster_prop, fill = species)) +
  geom_jitter(aes(col = species), width = 0.1, height = 0.1) +
  scale_color_manual(values = c(alpha("#ff9ebb", 0.7), alpha("#b9375e", 0.7))) +
  scale_fill_manual(values = c(alpha("#ff9ebb", 0.7), alpha("#b9375e", 0.7))) +
  
  ggnewscale::new_scale_color() +
  
  geom_smooth(aes(col = species), method = "lm") +
  
  scale_color_manual(values = c("#ff9ebb", "#b9375e"), guide = "none") +
  scale_fill_manual(values = c("#ff9ebb", "#b9375e"), guide = "none") +
  
  # scale_y_log10() +
  
  labs(y = "Proportion of cluster", x = "Number of genes") +
  
  theme_bw()


####################################################
#     ANALYSE DATA WITH SAMap CONNECTING PAIRS     #
####################################################

# # compute longer samap table
# ScilSlac_samapGenePairs_longer <- ScilSlac_samapGenePairs %>%
#   filter(pval1 < 0.05 & pval2 < 0.05) %>%
#   select(-c(pval1, pval2)) %>%
#   separate_rows(gene_pair, sep = ";") %>%
#   separate(mapping_clusters, sep = "\\.", remove = FALSE,
#            into = c("Scil_cluster", "Slac_cluster")) %>%
#   rename("mapping_id" = "mapping_clusters",
#          "gene_id" = "gene_pair") %>%
#   mutate(gene_id = str_remove(gene_id, "^.+_"))
# ScilSlac_samapGenePairs_longer

# calculate the proportion of connecting pairs in each og_type
contingency_connectingGenes <- orthogroup_full_set_with_geneCounts_clusterCounts_connections %>%
  distinct(orthogroup_id, og_type, is_connector) %>%
  count(og_type, is_connector) %>%
  complete(og_type, is_connector, fill = list(n = 0))
contingency_connectingGenes

# plot contingency as a matrix
contingency_connectingGenes_tile <- contingency_connectingGenes %>%
  ggplot(aes(og_type, is_connector, fill = n)) +
  geom_tile() +
  geom_text(aes(label = n), col = "white", fontface = "bold") +
  coord_fixed(ratio = 1)
contingency_connectingGenes_tile

ggsave("17_comparisons_1to1_orthologues/contingency_connectingGenes_tile.png",
       contingency_connectingGenes_tile, device = "png",
       height = 6, width = 10, dpi = 300, unit = "in", bg = "white")
ggsave("17_comparisons_1to1_orthologues/contingency_connectingGenes_tile.pdf",
       contingency_connectingGenes_tile, device = cairo_pdf,
       height = 6, width = 10, dpi = 300, unit = "in", bg = "white")

# plot contingency as barplots
contingency_connectingGenes_barplot <- contingency_connectingGenes %>%
  ggplot(aes(og_type, n, fill = is_connector)) +
  geom_col() +
  theme_for_plots
contingency_connectingGenes_barplot

ggsave("17_comparisons_1to1_orthologues/contingency_connectingGenes_barplot.png",
       contingency_connectingGenes_barplot, device = "png",
       height = 6, width = 10, dpi = 300, unit = "in", bg = "white")
ggsave("17_comparisons_1to1_orthologues/contingency_connectingGenes_barplot.pdf",
       contingency_connectingGenes_barplot, device = cairo_pdf,
       height = 6, width = 10, dpi = 300, unit = "in", bg = "white")

# compute fisher's exact
contingency_connectors %>%
  pivot_wider(names_from = is_connector, values_from = n, values_fill = 0) %>%
  column_to_rownames("og_type") %>%
  as.matrix() %>%
  fisher.test(simulate.p.value = TRUE, B = 10000)

# compute pairwisefisher's
contingency_connectors %>%
  pivot_wider(names_from = is_connector, values_from = n, values_fill = 0) %>%
  column_to_rownames("og_type") %>%
  as.matrix() %>%
  rstatix::pairwise_fisher_test(p.adjust.method = "BH")

# lookup table for OG names and types
og_type_lookup <- orthogroup_full_set_with_geneCounts_clusterCounts %>%
  distinct(orthogroup_id, og_type)
boxplots_per_og_type

# calculate for each leiden cluster the proportion of og_types
cluster_composition <- orthogroup_full_set_withClusterExpression %>%
  filter(!is.na(leiden_clusters)) %>% 
  distinct(orthogroup_id, leiden_clusters) %>%
  left_join(og_type_lookup, by = "orthogroup_id") %>%
  count(leiden_clusters, og_type) %>%
  group_by(leiden_clusters) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()
cluster_composition

# barplot of proportions of OGs type for each leiden cluster
cluster_composition %>%
  ggplot(aes(x = leiden_clusters, y = prop, fill = og_type)) +
  geom_col() +
  labs(y = "Proportion of orthogroups", x = "Leiden cluster") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# compute the number of exclusive 1 to 1 OGs for each cluster
cluster_specific_1to1 <- orthogroup_full_set_withClusterExpression %>%
  filter(!is.na(leiden_clusters)) %>%
  left_join(og_type_lookup, by = "orthogroup_id") %>%
  
  # calculate the total number of markers per cluster
  group_by(leiden_clusters) %>%
  mutate(total_markers_perCluster = n_distinct(orthogroup_id)) %>%
  ungroup() %>%
  
  filter(og_type == "1 to 1") %>%
  
  # keep only orthogroup/species combos that mark exactly 1 cluster
  group_by(orthogroup_id, species) %>%
  filter(n_distinct(leiden_clusters) == 1) %>%
  ungroup() %>%
  
  group_by(leiden_clusters) %>%
  mutate(n_1to1_specific_perCluster = n_distinct(orthogroup_id),
         prop_1to1_specific_perCluster = n_1to1_specific_perCluster/total_markers_perCluster) %>%
  ungroup() %>%
  distinct(leiden_clusters, prop_1to1_specific_perCluster, n_1to1_specific_perCluster, total_markers_perCluster)
cluster_specific_1to1

# subset the integrated object to only scil
scil_samap <- ScilSlac_samap %>%
  subset(species == "Scil")

# update metadata with the proportion of exclusive 1-to-1 OGs
scil_samap[[]] <- scil_samap[[]] %>%
  rownames_to_column() %>%
  mutate(Scil_leiden_clusters = paste0("Scil_", Scil_leiden_clusters)) %>%
  left_join(cluster_specific_1to1,
            by = join_by("Scil_leiden_clusters" == "leiden_clusters")) %>%
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
            by = join_by("Slac_leiden_clusters" == "leiden_clusters")) %>%
  replace(is.na(.), 0) %>%
  column_to_rownames()
slac_samap[[]] %>% View


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
  slice_sample(prop = 1) %>%
  mutate(species = case_when(species == "Scil" ~ "S. ciliatum",
                             species == "Slac" ~ "S. lacustris",
                             TRUE ~ species)) %>%
  
  ggplot(aes(UMAP_1, UMAP_2, col = species)) +
  geom_point(size = 0.8) +
  
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


scil_umap_prop_1to1 <- scil_umap_prop_1to1_raw@data %>%
  arrange(prop_1to1_specific_perCluster) %>%
  
  ggplot(aes(UMAP_1, UMAP_2, col = prop_1to1_specific_perCluster)) +
  geom_point(size = 0.8) +
  
  labs(title = expression(paste(bolditalic("S. ciliatum"), bold(" cell subset"))),
       col = "Transcriptional uniqueness\n(1-to-1 orthologues)") +
  
  scale_color_distiller(palette = "Blues", direction = 1, breaks = c(0, 0.25, 0.5)) +
  guides(size = guide_legend(label.position = "bottom"),
         color = guide_colorbar(barheight = 1)) +
  
  theme_bw(base_size = 12) +
  theme_for_UMAPS +
  theme(legend.position = "bottom",
        legend.title = element_text(hjust = 0.5, size = 10),
        legend.title.position = "top")
scil_umap_prop_1to1


slac_umap_prop_1to1 <- slac_umap_prop_1to1_raw@data %>%
  arrange(prop_1to1_specific_perCluster) %>%
  
  ggplot(aes(UMAP_1, UMAP_2, col = prop_1to1_specific_perCluster)) +
  geom_point(size = 0.8) +
  
  labs(title = expression(paste(bolditalic("S. lacustris"), bold(" cell subset"))),
       col = "Transcriptional uniqueness\n(1-to-1 orthologues)") +
  
  scale_color_distiller(palette = "Blues", direction = 1, breaks = c(0, 0.03, 0.06)) +
  guides(size = guide_legend(label.position = "bottom"),
         color = guide_colorbar(barheight = 1)) +
  
  theme_bw(base_size = 12) +
  theme_for_UMAPS +
  theme(legend.position = "bottom",
        legend.title = element_text(hjust = 0.5, size = 10),
        legend.title.position = "top")
slac_umap_prop_1to1

panel_umap <- scilslac_umap_integrated + scil_umap_prop_1to1 + slac_umap_prop_1to1
panel_umap

ggsave("17_comparisons_1to1_orthologues/panel_cluster_uniqueness.png",
       panel_umap, device = "png",
       height = 4*1.5, width = 10*1.5, dpi = 300, unit = "in", bg = "white")
ggsave("17_comparisons_1to1_orthologues/panel_cluster_uniqueness.pdf",
       panel_umap, device = cairo_pdf,
       height = 4*1.5, width = 10*1.5, dpi = 300, unit = "in", bg = "white")


cluster_specific_1toMany <- orthogroup_full_set_withClusterExpression %>%
  filter(!is.na(leiden_clusters)) %>%
  left_join(og_type_lookup, by = "orthogroup_id") %>%
  
  # calculate the total number of markers per cluster
  group_by(leiden_clusters) %>%
  mutate(total_markers_perCluster = n_distinct(orthogroup_id)) %>%
  ungroup() %>%
  
  filter(og_type == "Single copy\nin Sycon") %>%
  
  # keep only orthogroup/species combos that mark exactly 1 cluster
  group_by(orthogroup_id, species) %>%
  filter(n_distinct(leiden_clusters) == 1) %>%
  ungroup() %>%
  
  group_by(leiden_clusters) %>%
  mutate(n_1toMany_specific_perCluster = n_distinct(orthogroup_id),
         prop_1toMany_specific_perCluster = n_1toMany_specific_perCluster/total_markers_perCluster) %>%
  ungroup() %>%
  distinct(leiden_clusters, prop_1toMany_specific_perCluster, n_1toMany_specific_perCluster, total_markers_perCluster)
cluster_specific_1toMany

scil_samap[[]] <- scil_samap[[]] %>%
  rownames_to_column() %>%
  left_join(cluster_specific_1toMany,
            by = join_by("Scil_leiden_clusters" == "leiden_clusters")) %>%
  replace(is.na(.), 0) %>%
  column_to_rownames()

scil_samap %>%
  DimPlot2(features = "prop_1toMany_specific_perCluster", pt.size = 0.8,
           theme = theme_umap_arrows())

ScilSlac_samap %>%
  subset(species == "Scil") %>%
  DimPlot2(group.by = "Scil_leiden_clusters", label = TRUE, repel = TRUE)

scil_samap[[]] %>%
  select(Scil_leiden_clusters, Scil_seurat_clusters) %>%
  left_join(read.table("07_notableGenes_clusterAnnotation/cluster_identity.tsv", header = TRUE) %>%
              mutate(cluster_ID = as.factor(cluster_ID)),
            by = join_by("Scil_seurat_clusters" == "cluster_ID")) %>%
  count(Scil_leiden_clusters, Scil_seurat_clusters, value) %>%
  arrange(Scil_leiden_clusters, desc(n)) %>% View

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
    by = join_by(Slac_seurat_clusters_1 == rowname)
  ) %>%
  select(-Slac_seurat_clusters_2) %>% 
  count(Slac_leiden_clusters, Slac_seurat_clusters_1, cell_type) %>%
  arrange(Slac_leiden_clusters, desc(n)) %>%
  filter(Slac_leiden_clusters %in% c("Slac_26", "Slac_21", "Slac_6", "Slac_20", "Slac_0", "Slac_23", "Slac_25")) %>%
  View