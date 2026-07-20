#!/usr/bin/env Rscript

setwd("/data/evassvis/fn76/sycon/sycon_clusterAnnotation/sycon_cluster_analyses/")

library(tidyverse)
library(Seurat)
library(SeuratExtend)


pca_table <- read.table("17_comparisons_1to1_orthologues/orthogroup_analysis.tsv",
                        header = TRUE, sep = "\t", na.strings = c("", ".")) %>%
  as_tibble()


#############################################################
#     BOXPLOTS WITH PCA ORTHOGROUPS WITH EITHER SPECIES     #
#############################################################

# TOTAL COUNT: 5,460

df_clusters_per_genes <- pca_table %>%
  mutate(unique_ID = paste0(orthogroup, "_", pca_id)) %>% 
  
  filter(pca_kind %in% c("minimal", "maximal+minimal"),
         !(is.na(Sycon_ids) & is.na(Spongilla_ids))) %>%
  
  select(unique_ID,
         Sycon_id_count, Sycon_clusters,
         Spongilla_id_count, Spongilla_clusters) %>%
  
  mutate(across(everything(), ~replace_na(., "0")),
         across(everything(), ~as.character(.))) %>%
  
  pivot_longer(-unique_ID) %>%
  
  separate(name, into = c("species", "group"), sep = "_", extra = "merge") %>%
  
  mutate(group = str_replace(group, "id_count", "number_genes"),
         count = case_when(
           group == "clusters" & value != "0" ~ str_count(value, ";") + 1,
           group == "clusters" & value == "0" ~ 0,
           group == "number_genes" ~ as.numeric(value)
           )) %>%
  
  select(-value) %>%
  
  pivot_wider(names_from = group, values_from = count) %>%
  
  filter(clusters > 0)
  
  
boxplot_clusters_per_genes <- df_clusters_per_genes %>%
  ggplot(aes(as_factor(clusters), number_genes, fill = species)) +
  geom_point(aes(color = species),
             position = position_jitterdodge(jitter.width = 0.3,
                                             jitter.height = 0.3,
                                             dodge.width = 1)) +
  geom_boxplot(outliers = FALSE, 
               position = position_dodge(width = 1)) +
  
  ggpubr::stat_compare_means(aes(group = species), label = "p.signif") +
  
  labs(x = "Number of cell clusters", y = "Number of genes") +
  
  theme_bw()
boxplot_clusters_per_genes

ggsave("17_comparisons_1to1_orthologues/boxplot_clusters_per_genes.png",
       boxplot_clusters_per_genes, device = "png",
       width = 10, height = 6, dpi = 300, unit = "in", bg = "white")
ggsave("17_comparisons_1to1_orthologues/boxplot_clusters_per_genes.pdf",
       boxplot_clusters_per_genes, device = cairo_pdf,
       width = 10, height = 6, dpi = 300, unit = "in", bg = "white")


#######################################################################
#     BOXPLOTS WITH PCA ORTHOGROUPS WITH BOTH SINGLE-COPY SPECIES     #
#######################################################################

sycon_singleCopy_list <- df_clusters_per_genes %>% 
  filter(species == "Sycon" & number_genes == "1") %>%
  pull(unique_ID)

length(sycon_singleCopy_list)

df_clusters_per_genes %>% 
  filter(unique_ID %in% sycon_singleCopy_list) %>% 
  
  ggplot(aes(x = species, y = clusters, group = unique_ID)) +
  
  # geom_line(alpha = 0.3, color = "grey50") +
  
  geom_point(aes(color = species), size = 2) +
  
  ggpubr::stat_compare_means(method = "wilcox.test", paired = TRUE) +
  
  scale_y_continuous(limits = c(0, 5), breaks = seq(0, 5)) +
  
  labs(x = "Species", y = "Number of clusters") +
  
  theme_bw()

spongilla_singleCopy_list <- df_clusters_per_genes %>% 
  filter(species == "Spongilla" & number_genes == "1") %>%
  pull(unique_ID)

df_clusters_per_genes %>% 
  filter(unique_ID %in% spongilla_singleCopy_list) %>% 
  
  ggplot(aes(as_factor(species), clusters, fill = species)) +
  geom_point(aes(color = species),
             position = position_jitterdodge(dodge.width = 1)) +
  geom_boxplot(outliers = FALSE, 
               position = position_dodge(width = 1)) +
  
  ggpubr::stat_compare_means(label = "p.signif") +
  
  scale_y_continuous(limits = c(0, 5), breaks = seq(0, 5)) +
  
  labs(x = "Species", y = "Number of clusters") +
  
  theme_bw()

single_copy_comparison <- df_clusters_per_genes %>% 
  complete(unique_ID, species, explicit = FALSE, fill = list(number_genes = 0, clusters = 0)) %>%
  group_by(unique_ID) %>%
  mutate(category = case_when(all(number_genes == 1) ~ "Single-copy OGs\nin both species",
                              any(species == "Sycon" & number_genes == 1) ~ "Single-copy OGs\nin Sycon",
                              any(species == "Spongilla" & number_genes == 1) ~ "Single-copy OGs\nin Spongilla",
                              TRUE ~ "Neither")) %>%
  ungroup() %>%

  filter(category != "Neither") %>%
  
  ggplot(aes(x = species,
             y = clusters)) +
  
  # geom_line(alpha = 0.3) +
  geom_jitter(width = 0.1, height = 0.1) +
  # geom_violin(aes(fill = category)) +
  
  guides(fill = "none") +
  
  facet_wrap(~category) +
  
  # ggpubr::stat_compare_means(method = "wilcox.test", paired = TRUE) +
  
  scale_y_continuous(limits = c(0, 5), breaks = seq(0, 5)) +
  
  theme_bw()
single_copy_comparison

ggsave("17_comparisons_1to1_orthologues/single_copy_comparison.png",
       single_copy_comparison, device = "png",
       width = 10, height = 6, dpi = 300, unit = "in", bg = "white")
ggsave("17_comparisons_1to1_orthologues/single_copy_comparison.pdf",
       single_copy_comparison, device = cairo_pdf,
       width = 10, height = 6, dpi = 300, unit = "in", bg = "white")


###########################################
#     ANALYSIS WITHIN THE SAMap SPACE     #
###########################################

# before starting, load slac gene name conversion table
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


ScilSlac_samap <- sceasy::convertFormat("05_SAMap_porifera/ScilSlac_leiden3Clusters_samap.h5ad",
                                        from = "anndata", to = "seurat")

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



scil_samap <- ScilSlac_samap %>%
  subset(species == "Scil") %>%
  SCTransform()
slac_samap <- ScilSlac_samap %>%
  subset(species == "Slac") %>%
  SCTransform()

scil_samap_allMarkers <- scil_samap %>%
  FindAllMarkers(group.by = "Scil_leiden_clusters") %>%
  rename("leiden_clusters" = "cluster") %>%
  mutate(gene = str_replace(gene, "Scil-", ""))

write.table(scil_samap_allMarkers, "17_comparisons_1to1_orthologues/scil_samap_allMarkers_leidenClusters.tsv",
            col.names = TRUE, sep = "\t", quote = FALSE, row.names = FALSE)

slac_samap_allMarkers <- slac_samap %>%
  FindAllMarkers(group.by = "Slac_leiden_clusters") %>%
  rename("leiden_clusters" = "cluster") %>%
  mutate(gene = str_replace(gene, "Slac-", "")) %>%
  left_join(slac_gene_names_conversion, by = join_by("gene" == "gene_id"),
            relationship = "many-to-many")

write.table(slac_samap_allMarkers, "17_comparisons_1to1_orthologues/slac_samap_allMarkers_leidenClusters.tsv",
            col.names = TRUE, sep = "\t", quote = FALSE, row.names = FALSE)

df_leidenClusters_per_genes <- pca_table %>%
  mutate(unique_ID = paste0(orthogroup, "_", pca_id)) %>% 
  
  filter(pca_kind %in% c("minimal", "maximal+minimal")) %>% 
  
  select(unique_ID,
         Sycon_ids, Spongilla_ids,
         Sycon_id_count, Spongilla_id_count) %>%
  
  pivot_longer(c(Sycon_ids, Spongilla_ids),
               names_to = "species", values_to = "genes") %>%
  
  mutate(species = str_replace(species, "_ids", "")) %>%
  
  separate_rows(genes, sep = ";") %>%
  
  mutate(genes = str_replace(genes, "\\.[0-9]+$", "")) %>%
  
  left_join(scil_samap_allMarkers %>%
              add_row(slac_samap_allMarkers %>%
                        select(-c(gene, transcript_id)) %>%
                        rename("gene" = "protein_id")) %>%
              select(leiden_clusters, gene),
            by = join_by("genes" == "gene"),
            relationship = "many-to-many") %>%

  group_by(unique_ID, species) %>%
  # drop_na() %>% # this would remove rows with 0 count on genes and on leiden_clusters
  summarise(n_genes = n_distinct(genes[!is.na(genes)]),
            n_leiden_clusters = n_distinct(leiden_clusters[!is.na(leiden_clusters)]),
            .groups = "drop")
 
df_leidenClusters_per_genes_withCategories <- df_leidenClusters_per_genes %>% 
  complete(unique_ID, species, explicit = FALSE, fill = list(n_genes = 0, n_leiden_clusters = 0)) %>%
  group_by(unique_ID) %>% 
  mutate(category = case_when(all(n_genes == 1) ~ "1-to-1",
                              any(species == "Sycon" & n_genes == 0) ~ "Species specific",
                              any(species == "Spongilla" & n_genes == 0) ~ "Species specific",
                              any(species == "Sycon" & n_genes == 1) ~ "many-to-1\n(Spongilla-Sycon)",
                              any(species == "Spongilla" & n_genes == 1) ~ "1-to-many\n(Spongilla-Sycon)",
                              TRUE ~ "many-to-many")) %>%
  ungroup()



leidenClusters_per_genes <- df_leidenClusters_per_genes_withCategories %>%
  # filter() %>%
  # mutate(n_genes = n_genes + 0.00000001) %>%
  ggplot(aes(n_leiden_clusters, n_genes, fill = species)) +
  geom_jitter(aes(col = species), width = 0.3, height = 0.3) +
  scale_color_manual(values = c("#9aabe6", "#e3b498")) +
  
  ggnewscale::new_scale_color() +
  
  geom_smooth(aes(col = species), method = "lm") +
  
  # geom_point(aes(color = species),
  #            position = position_jitterdodge(jitter.width = 0.3,
  #                                            jitter.height = 0.3,
  #                                            dodge.width = 1)) +
  # geom_boxplot(outliers = FALSE, 
  #              position = position_dodge(width = 1)) +
  
  scale_color_manual(values = c("#4a6ee2", "#dc6c29")) +
  scale_fill_manual(values = c("#9aabe6", "#e3b498")) +
  
  scale_y_log10() +
  
  labs(x = "Number of SAMap leiden clusters", y = "Number of genes") +
  
  theme_bw()
leidenClusters_per_genes

model <- lm(n_genes ~ n_leiden_clusters * species,
            data = df_leidenClusters_per_genes_withCategories)
summary(model)
broom::tidy(model)

ggsave("17_comparisons_1to1_orthologues/leidenClusters_per_genes_fit.png",
       leidenClusters_per_genes, device = "png",
       width = 7, height = 4, dpi = 300, unit = "in", bg = "white")
ggsave("17_comparisons_1to1_orthologues/leidenClusters_per_genes_fit.pdf",
       leidenClusters_per_genes, device = cairo_pdf,
       width = 7, height = 4, dpi = 300, unit = "in", bg = "white")

df_leidenClusters_per_genes %>%
  ggplot(aes(as_factor(n_genes), n_leiden_clusters, fill = species)) +
  geom_point(aes(color = species),
             position = position_jitterdodge(jitter.width = 0.3,
                                             jitter.height = 0.3,
                                             dodge.width = 0.6)) +  
  geom_boxplot(outliers = FALSE, width = 0.5,
               position = position_dodge(width = 0.6)) +
  
  ggpubr::stat_compare_means(aes(group = species), label = "p.signif") +
  
  # labs(x = "Number of SAMap leiden clusters", y = "Number of genes") +
  
  theme_bw()

boxplot_orthogroup_categories <- df_leidenClusters_per_genes_withCategories %>%
  
  filter(n_leiden_clusters > 0) %>%
  
  ggplot(aes(x = category, y = n_leiden_clusters,
             fill = species)) +
  
  # geom_violin(aes(fill = species)) +
  geom_point(aes(color = species),
             position = position_jitterdodge(jitter.width = 0.3,
                                             jitter.height = 0.3,
                                             dodge.width = 0.6)) +  
  geom_boxplot(outliers = FALSE, width = 0.5,
               position = position_dodge(width = 0.6)) +
  
  scale_fill_manual(values = c("#4a6ee2", "#dc6c29")) +
  scale_color_manual(values = c("#9aabe6", "#e3b498")) +
  
  labs(x = "Orthogroup category", y = "# leiden clusters") +
  
  ggpubr::stat_compare_means(aes(group = species), label = "p.signif") +
  
  theme_bw()

boxplot_orthogroup_categories


ggsave("17_comparisons_1to1_orthologues/boxplot_orthogroup_categories.png",
       boxplot_orthogroup_categories, device = "png",
       height = 4*1.5, width = 7*1.5, dpi = 300, unit = "in", bg = "white")
ggsave("17_comparisons_1to1_orthologues/boxplot_orthogroup_categories.pdf",
       boxplot_orthogroup_categories, device = cairo_pdf,
       height = 4*1.5, width = 7*1.5, dpi = 300, unit = "in", bg = "white")


#################
#     UCell     #
#################

oneToMany_spongilla_sycon <- df_leidenClusters_per_genes_withCategories %>%
  complete(unique_ID, species, explicit = FALSE, fill = list(n_genes = 0, n_leiden_clusters = 0)) %>%
  group_by(unique_ID) %>%
  mutate(category = case_when(all(n_genes == 1) ~ "1-to-1",
                              any(species == "Sycon" & n_genes == 0) ~ "Spongilla specific",
                              any(species == "Spongilla" & n_genes == 0) ~ "Sycon specific",
                              any(species == "Sycon" & n_genes == 1) ~ "many-to-1\n(Spongilla-Sycon)",
                              any(species == "Spongilla" & n_genes == 1) ~ "1-to-many\n(Spongilla-Sycon)",
                              TRUE ~ "many-to-many")) %>%
  filter(category == "1-to-many\n(Spongilla-Sycon)",
         species == "Sycon",
         n_leiden_clusters == 1) %>%
  arrange(desc(n_leiden_clusters))

manyToOne_spongilla_sycon <- df_leidenClusters_per_genes_withCategories %>%
  complete(unique_ID, species, explicit = FALSE, fill = list(n_genes = 0, n_leiden_clusters = 0)) %>%
  group_by(unique_ID) %>%
  mutate(category = case_when(all(n_genes == 1) ~ "1-to-1",
                              any(species == "Sycon" & n_genes == 0) ~ "Spongilla specific",
                              any(species == "Spongilla" & n_genes == 0) ~ "Sycon specific",
                              any(species == "Sycon" & n_genes == 1) ~ "many-to-1\n(Spongilla-Sycon)",
                              any(species == "Spongilla" & n_genes == 1) ~ "1-to-many\n(Spongilla-Sycon)",
                              TRUE ~ "many-to-many")) %>%
  filter(category == "many-to-1\n(Spongilla-Sycon)",
         species == "Spongilla",
         n_leiden_clusters == 1) %>%
  arrange(desc(n_leiden_clusters))

oneToMany_spongilla_sycon %>% arrange(desc(n_genes))
# OG0000708_PCA001 (n_gene = 5; n_leiden_clusters = 26)

manyToOne_spongilla_sycon %>% arrange(desc(n_genes)) %>% View
# OG0000224_PCA002 (n_genes = 9; n_leiden_clusters = 30)

spongilla_genes <- pca_table %>%
  filter(orthogroup == "OG0001045" & pca_id == "PCA001") %>%
  pull(Spongilla_ids) %>%
  strsplit(., ";") %>%
  .[[1]] %>%
  paste0("Slac-", .) %>%
  list()

sycon_genes <- pca_table %>%
  filter(orthogroup == "OG0001045" & pca_id == "PCA001") %>%
  pull(Sycon_ids) %>%
  strsplit(., ";") %>%
  .[[1]] %>%
  paste0("Scil-", .) %>%
  list()
scil_samap_allMarkers %>% filter(gene == "g9334")

sycon_UCell <- AddModuleScore_UCell(ScilSlac_samap,
                                    features = sycon_genes)
spongilla_UCell <- AddModuleScore_UCell(slac_samap,
                                        features = spongilla_genes)

sycon_UCell_plot <- SeuratExtend::DimPlot2(sycon_UCell, features = "signature_1_UCell")
spongilla_UCell_plot <- SeuratExtend::DimPlot2(spongilla_UCell[[]], features = "signature_1_UCell")

