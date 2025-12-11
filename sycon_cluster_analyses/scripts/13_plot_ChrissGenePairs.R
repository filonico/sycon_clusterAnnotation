setwd("/lustre/alice3/data/evassvis/fn76/sycon/sycon_clusterAnnotation/sycon_cluster_analyses")

library(Seurat)
library(SeuratExtend)
library(reticulate)
library(scCustomize)
library(tidyverse)
library(metacell)


#######################
#     LOAD INPUTS     #
#######################

spongilla_genename_conversion <- read.table("00_input/slac_genome_genename_conversion.ls", sep = "_")

# import Chris's gene pairs
OneToOne_orthologs <- read.table("00_input/in_paralogs.tsv", header = FALSE, sep = "\t", na.strings = "") %>%
  drop_na() %>%
  select(V8, V9) %>%
  filter(!grepl(";", V8),
         !grepl(";", V9)) %>%
  mutate(V9 = str_replace_all(V9, "\\.[0-9]$", "")) %>%
  left_join(spongilla_genename_conversion, by = join_by("V9" == "V1")) %>%
  mutate(gene_pairs = paste0("Scil_", V8, ";Slac_", V2),
         OneToOne = "yes")
OneToOne_orthologs

# read the SAMap identified gene pairs
samap_gene_pairs <- read.table("05_SAMap_porifera/02_gene_pairs/ScilSlac_leiden3Clusters_all_samapGenePairs.tsv",
           header = TRUE, sep = "\t", na.strings = "") %>%
  select(where(is.character)) %>%
  pivot_longer(everything(),
               names_to = "cell_pairs",
               values_to = "gene_pairs") %>%
  drop_na()

barplot <- samap_gene_pairs %>%
  left_join(OneToOne_orthologs) %>%
  select(-c(V8, V9, V2)) %>%
  replace_na(list(OneToOne = "No")) %>%
  
  ggplot(aes(x = cell_pairs, fill = OneToOne)) +
  geom_bar() +
  
  scale_fill_manual(values = c("#e69f00", "#2271b2")) +
  
  labs(x = "Mapped cell-cluster pairs", y = "Number of genes connecting the pairs") +
  
  theme_bw(base_size = 12) %+replace%
  theme(plot.background = element_rect(fill = "transparent", colour = NA), 
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey90", lineend = "round"),
        panel.border = element_blank(),
        legend.background = element_rect(fill = "transparent", colour = NA),
        legend.key = element_rect(fill = "transparent", colour = NA),
        legend.key.width = unit(0.4, "cm"),
        legend.key.height = unit(0.4, "cm"),
        legend.position = "bottom",
        plot.title = element_text(size = 13, hjust = 0.0, vjust = 1.75, face = "bold"),
        axis.line = element_line(color = "black", linewidth = 0.6, lineend = "round"),
        axis.ticks = element_line(colour = "black", linewidth = 0.6, lineend = "round"),
        axis.ticks.length = unit(0.20, "cm"),
        axis.text.y = element_text(color = "black", hjust = 1,
                                   margin = margin(t = 0, r = 4, b = 0, l = 10)),
        axis.text.x = element_text(color = "black", angle = 45, hjust = 1, vjust = 1,
                                   margin = margin(t = 4, r = 0, b = 10, l = 0)),
        axis.title.x = element_text(size = 13, angle = 0,
                                    margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.y = element_text(size = 13, angle = 90,
                                    margin = margin(t = 10, r = 0, b = 0, l = 0)))
barplot

ggsave("07_notableGenes_clusterAnnotation/SAMap_genPairs_ChrissOneToOne.pdf",
       barplot, device = cairo_pdf,
       height = 6, width = 8, dpi = 300, bg = "white")
ggsave("07_notableGenes_clusterAnnotation/SAMap_genPairs_ChrissOneToOne.png",
       barplot, device = "png",
       height = 6, width = 8, dpi = 300, bg = "white")




# add a meta.data column with merged cell annotation
ScilSlac_samap@meta.data <- ScilSlac_samap@meta.data %>%
  mutate(merged = pmap_chr(list(Slac_cell_type, Scil_seurat_clusters), ~ {
    if (..1 != "unassigned") {
      paste0("Slac_", ..1)
    } else if (..2 != "unassigned") {
      paste0("Scil_", ..2)
    } else {
      NA_character_
    }
  }),
  species_full = str_replace(species_full, "Spongilla lacustris", "Spongilla"),
  species_full = str_replace(species_full, "Sycon ciliatum", "Sycon"))

# add a meta.data column with merged cell annotation
AqueSlac_samap@meta.data <- AqueSlac_samap@meta.data %>%
  mutate(merged = pmap_chr(list(Slac_cell_type, Aque_cell_type), ~ {
    if (..1 != "unassigned") {
      paste0("Slac_", ..1)
    } else if (..2 != "unassigned") {
      paste0("Aque_", ..2)
    } else {
      NA_character_
    }
  }),
  species_full = str_replace(species_full, "Spongilla lacustris", "Spongilla"),
  species_full = str_replace(species_full, "Amphimedon queenslandica", "Amphimedon"))

# add a meta.data column with merged cell annotation
AqueScil_samap@meta.data <- AqueScil_samap@meta.data %>%
  mutate(merged = pmap_chr(list(Scil_seurat_clusters, Aque_cell_type), ~ {
    if (..1 != "unassigned") {
      paste0("Scil_", ..1)
    } else if (..2 != "unassigned") {
      paste0("Aque_", ..2)
    } else {
      NA_character_
    }
  }),
  species_full = str_replace(species_full, "Sycon ciliatum", "Sycon"),
  species_full = str_replace(species_full, "Amphimedon queenslandica", "Amphimedon"))


# plot the UMAP with species assignment
ScilSlac_umap_species <- ScilSlac_samap %>%
  SetIdent(value = 'species_full') %>%
  DimPlot2(theme = theme_umap_arrows(), pt.size = 0.8, cols = c("#ffa500", "#ff186e")) +
  theme(legend.position = "inside",
        legend.position.inside = c(0,0))

ScilSlac_umap_species

AqueSlac_umap_species <- AqueSlac_samap %>%
  SetIdent(value = 'species_full') %>%
  DimPlot(theme = theme_umap_arrows(), alpha = 0.5, pt.size = 2, cols = c("#2690f8", "#ffa500")) +
  theme(legend.position = "inside",
        legend.position.inside = c(0,0))

AqueSlac_umap_species

AqueScil_umap_species <- AqueScil_samap %>%
  SetIdent(value = 'species_full') %>%
  DimPlot2(theme = theme_umap_arrows(), pt.size = 2, cols = c("#2690f8", "#ff186e")) +
  theme(legend.position = "inside",
        legend.position.inside = c(0,0))

AqueScil_umap_species

# plot the UMAP with cell-type assignment
ScilSlac_umap_cellAnno <- DimPlot(object = ScilSlac_samap, group.by = 'merged', alpha = 0.5, pt.size = 0.8) %>%
  LabelClusters(id = 'merged', repel = TRUE, max.overlaps = Inf, size = 3) +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())
ScilSlac_umap_cellAnno

AqueSlac_umap_cellAnno <- DimPlot(object = AqueSlac_samap, group.by = 'merged', alpha = 0.5, pt.size = 2) %>%
  LabelClusters(id = 'merged', repel = TRUE, max.overlaps = Inf, size = 3) +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())
AqueSlac_umap_cellAnno

AqueScil_umap_cellAnno <- DimPlot(object = AqueScil_samap, group.by = 'merged', alpha = 0.5, pt.size = 2) %>%
  LabelClusters(id = 'merged', repel = TRUE, max.overlaps = Inf, size = 3) +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())
AqueScil_umap_cellAnno

# identify common features between Chris's gene pairs and SAMap analysed data
common_features <- as.vector(t(ScilSlac_testPairs %>% select(Slac,Scil))) %in% rownames(ScilSlac_samap)

# identify features not present in the SAMap analysed data
not_present <- as.vector(t(ScilSlac_testPairs %>% select(Slac,Scil)))[!common_features]

# downsample Chris's gene pairs to only those present in the SAMap analysis
ScilSlac_testPairs <- ScilSlac_testPairs %>%
  filter(!Slac %in% not_present)

# plot pairs of features from Scil and Slac and save it
ScilSlac_genePairs_UMAP <- ScilSlac_samap %>%
  FeaturePlot3.grid(color = "ryb", features = as.vector(t(ScilSlac_testPairs)),
                    pt.size = 0.8)
ggsave("05_SAMap/02_gene_pairs/plot_ChrissGenePairs.png",
       ScilSlac_genePairs_UMAP, device = png,
       dpi = 300, height = 14*3, width = 14*3, units = ("in"), bg = 'white')

# read the SAMap identified gene pairs
samap_gene_pairs <- read.table("05_SAMap/02_gene_pairs/ScilSlac_leiden3Clusters_all_samapGenePairs.tsv", header = TRUE, sep = "\t") %>%
  # let the gene ID match those from Chris's test pairs
  mutate(across(where(is.character), ~ str_replace_all(., "_", "-") %>%
                  str_remove_all("-g[0-9]$")))

# pivot longer samap gene pairs
samap_gene_pairs_longer <- samap_gene_pairs %>%
  select(where(is.character)) %>%
  pivot_longer(everything(),
               names_to = "cell_pairs",
               values_to = "gene_pairs") %>%
  mutate(across(where(is.character),  ~na_if(.x, ""))) %>%
  mutate(cell_pairs = str_replace(cell_pairs, "\\.", "|"),
         gene_pairs = str_replace(gene_pairs, ";", "|")) %>%
  drop_na(gene_pairs)

write.table(samap_gene_pairs_longer, file = "05_SAMap/02_gene_pairs/ScilSlac_leiden3Clusters_all_samapGenePairs_longer.tsv",
            quote = FALSE, sep = "\t", row.names = FALSE)

# count the number of gene pairs supporting cell type mappings
samap_cellPairs_count <- samap_gene_pairs_longer %>%
  count(cell_pairs) %>%
  arrange(desc(n))

write.table(samap_cellPairs_count, file = "05_SAMap/02_gene_pairs/ScilSlac_leiden3Clusters_all_samapGenePairs_counts.tsv",
            quote = FALSE, sep = "\t", row.names = FALSE)

# extract only the SAMap gene pairs matching those from Chris's data
connecting_gene_pairs <- samap_gene_pairs_longer %>%
  filter(gene_pairs %in% ScilSlac_testPairs$combined)

write.table(connecting_gene_pairs, file = "05_SAMap/02_gene_pairs/ScilSlac_leiden3Clusters_all_samapGenePairs_longer_chrissMatches.tsv",
            quote = FALSE, sep = "\t", row.names = FALSE)

connecting_gene_pairs_count <- connecting_gene_pairs %>%
  count(cell_pairs) %>%
  arrange(desc(n))

write.table(connecting_gene_pairs_count, file = "05_SAMap/02_gene_pairs/ScilSlac_leiden3Clusters_all_samapGenePairs_chrissMatches_counts.tsv",
            quote = FALSE, sep = "\t", row.names = FALSE)

df_to_plot <- merge(samap_cellPairs_count, connecting_gene_pairs_count, by = "cell_pairs", all = TRUE) %>%
  rename(c(total_pairs = n.x, matching_pairs = n.y)) %>%
  mutate(matching_pairs = replace_na(matching_pairs, 0),
         not_matching_pairs = total_pairs - matching_pairs) %>%
  select(-total_pairs) %>%
  pivot_longer(-cell_pairs, names_to = "pair_type", values_to = "count") %>%
  group_by(cell_pairs) %>%
  mutate(total = sum(count)) %>%
  ungroup() %>%
  mutate(cell_pairs = fct_reorder(cell_pairs, total, .desc = TRUE))

barplot_counts <- df_to_plot %>% 
  filter(pair_type != "total_pairs") %>%
  
  ggplot(aes(x = cell_pairs, y = count)) +
    geom_bar(aes(fill = pair_type), stat = "identity", width = 0.8) +
  
    colorspace::scale_fill_discrete_diverging(palette = "Cyan-Mage",
                                              labels = c("matching_pairs" = "Overlapping gene pairs",
                                                         "not_matching_pairs" = "SAMap-only gene pairs")) +
    labs(x = "Mapping cell types", y = "Counts", fill = "") +
  
    coord_cartesian(clip = "off") +
    scale_y_continuous(expand = c(0, 0), limits = c(-4, max(df_to_plot$count)+30)) +
    scale_x_discrete() +
  
    ggtitle("Number of gene pairs per mapping cell cluster") +
  
    geom_text(df_to_plot %>% distinct(cell_pairs, total),
              mapping = aes(x = cell_pairs, y = total, label = total), vjust = -0.5, size = 2) +
  
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text = element_text(size = 8),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.line = element_line(linewidth = 0.5),
          legend.background = element_rect(fill = 'transparent', color = NA),
          legend.position = "inside",
          legend.position.inside = c(0.8, 0.9),
          plot.margin = margin(.5, .5, .5, 2, "cm"))

barplot_counts

ggsave("05_SAMap/02_gene_pairs/ScilSlac_leiden3Clusters_all_samapGenePairs_barplot.png",
       barplot_counts, device = 'png',
       dpi = 300, height = 7, width = 7, units = ("in"), bg = 'white')

# connecting_gene_pairs %>%
#   filter(str_detect(gene_pairs, "8861"))

ScilSlac_samap@meta.data <- ScilSlac_samap@meta.data %>%
  mutate(integrated_cellAnnotation = merged) %>%
  mutate(integrated_cellAnnotation = str_replace_all(integrated_cellAnnotation, "Slac_Myopeptidocytes[0-9]$", "MYOPEPTIDOCYTES"),
         integrated_cellAnnotation = str_replace_all(integrated_cellAnnotation, "Scil_1$", "MYOPEPTIDOCYTES"),
         integrated_cellAnnotation = str_replace_all(integrated_cellAnnotation, "Slac_Mesocytes+$", "MESOCYTES"),
         integrated_cellAnnotation = str_replace_all(integrated_cellAnnotation, "Scil_2$", "MESOCYTES 1"),
         integrated_cellAnnotation = str_replace_all(integrated_cellAnnotation, "Scil_7$", "MESOCYTES 1"),
         integrated_cellAnnotation = str_replace_all(integrated_cellAnnotation, "Slac_Mesocytes 1$", "MESOCYTES 1"),
         integrated_cellAnnotation = str_replace_all(integrated_cellAnnotation, "Slac_Archaeocytes", "ARCHAEOCYTES"),
         integrated_cellAnnotation = str_replace_all(integrated_cellAnnotation, "Scil_4$", "ARCHAEOCYTES"),
         integrated_cellAnnotation = str_replace_all(integrated_cellAnnotation, "Scil_5$", "ARCHAEOCYTES"),
         integrated_cellAnnotation = str_replace_all(integrated_cellAnnotation, "Scil_10$", "ARCHAEOCYTES"),
         integrated_cellAnnotation = str_replace_all(integrated_cellAnnotation, "Scil_14$", "BASOPINACOCYTES"),
         integrated_cellAnnotation = str_replace_all(integrated_cellAnnotation, "Slac_Basopinacocytes$", "BASOPINACOCYTES"),
         integrated_cellAnnotation = str_replace_all(integrated_cellAnnotation, "Scil_15$", "METABOLOCYTES"),
         integrated_cellAnnotation = str_replace_all(integrated_cellAnnotation, "Slac_Metabolocytes[0-9]$", "METABOLOCYTES"),
         integrated_cellAnnotation = str_replace_all(integrated_cellAnnotation, "Slac_23$", "METABOLOCYTES"),
         integrated_cellAnnotation = str_replace_all(integrated_cellAnnotation, "Scil_19$", "SLAC6"),
         integrated_cellAnnotation = str_replace_all(integrated_cellAnnotation, "Slac_6$", "SLAC6"),
         integrated_cellAnnotation = str_replace_all(integrated_cellAnnotation, "Scil_21$", "SCLEROCYTES"),
         integrated_cellAnnotation = str_replace_all(integrated_cellAnnotation, "Scil_22$", "CHOANOBLASTS"),
         integrated_cellAnnotation = str_replace_all(integrated_cellAnnotation, "Slac_Choanoblasts[0-9]$", "CHOANOBLASTS"),
         integrated_cellAnnotation = str_replace_all(integrated_cellAnnotation, "Slac_Incurrent.Pinacocytes[0-9]$", "PINACOCYTES"),
         integrated_cellAnnotation = str_replace_all(integrated_cellAnnotation, "Scil_25$", "PINACOCYTES"),
         integrated_cellAnnotation = str_replace_all(integrated_cellAnnotation, "Scil_26$", "SCLEROCYTES"),
         integrated_cellAnnotation = str_replace_all(integrated_cellAnnotation, "Slac_Sclerocytes$", "SCLEROCYTES"),
         integrated_cellAnnotation = str_replace_all(integrated_cellAnnotation, "Scil_27$", "SLAC6"),
         integrated_cellAnnotation = str_replace_all(integrated_cellAnnotation, "Scil_28$", "MESOCYTES 3"),
         integrated_cellAnnotation = str_replace_all(integrated_cellAnnotation, "Slac_Mesocytes 3$", "MESOCYTES 3"),
         integrated_cellAnnotation = str_replace_all(integrated_cellAnnotation, "Slac_Apendopinacocytes[0-9]$", "Slac_Apendopinacocytes"),
         integrated_cellAnnotation = str_replace_all(integrated_cellAnnotation, "Scil_32$", "ARCHAEOCYTES")
  )

ScilSlac_samap <- ScilSlac_samap %>%
  SetIdent(value = "integrated_cellAnnotation")

ARCHAEOCYTES <- WhichCells(ScilSlac_samap, idents = "ARCHAEOCYTES")
BASOPINACOCYTES <- WhichCells(ScilSlac_samap, idents = "BASOPINACOCYTES")
CHOANOBLASTS <- WhichCells(ScilSlac_samap, idents = "CHOANOBLASTS")
MESOCYTES_1 <- WhichCells(ScilSlac_samap, idents = "MESOCYTES 1")
MESOCYTES_3 <- WhichCells(ScilSlac_samap, idents = "MESOCYTES 3")
METABOLOCYTES <- WhichCells(ScilSlac_samap, idents = "METABOLOCYTES")
MYOPEPTIDOCYTES <- WhichCells(ScilSlac_samap, idents = "MYOPEPTIDOCYTES")
PINACOCYTES <- WhichCells(ScilSlac_samap, idents = "PINACOCYTES")
SCLEROCYTES <- WhichCells(ScilSlac_samap, idents = "SCLEROCYTES")
SLAC6 <- WhichCells(ScilSlac_samap, idents = "SLAC6")

ScilSlac_umap_cellAnno_integrated <- DimPlot(object = ScilSlac_samap, group.by = 'integrated_cellAnnotation', alpha = 0.5, pt.size = 0.8, cols = "grey",
                                             cells.highlight = list(ARCHAEOCYTES,BASOPINACOCYTES,CHOANOBLASTS,MESOCYTES_3,MESOCYTES_1,METABOLOCYTES,MYOPEPTIDOCYTES,PINACOCYTES,SCLEROCYTES,SLAC6),
                                             cols.highlight = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF","#E69F00", "#00BA38", "#619CFF", "#F564E3", "#FF61CC", "#B79F00")) %>%
  LabelClusters(id = 'integrated_cellAnnotation', repel = TRUE, max.overlaps = Inf, size = 3) +
    theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())

ScilSlac_umap_cellAnno_integrated

# plot a panel with the 2 UMAPs and save it to a file
panel <- ggpubr::ggarrange(ScilSlac_umap_species, ScilSlac_umap_cellAnno, ScilSlac_umap_cellAnno_integrated, ncol = 3)
panel
panel <- ggpubr::ggarrange(ScilSlac_umap_species, ScilSlac_umap_cellAnno, ncol = 2)
panel
panel <- ggpubr::ggarrange(AqueSlac_umap_species, AqueSlac_umap_cellAnno, ncol = 2)
panel
panel <- ggpubr::ggarrange(AqueScil_umap_species, AqueScil_umap_cellAnno, ncol = 2)
panel

ggsave("panel_ScilSlac_samapIntegration.png",
       panel, device = png,
       dpi = 300, height = 8, width = 16, units = ("in"), bg = 'white')
ggsave("panel_AqueSlac_samapIntegration.png",
       panel, device = png,
       dpi = 300, height = 8, width = 16, units = ("in"), bg = 'white')
ggsave("panel_AqueScil_samapIntegration.png",
       panel, device = png,
       dpi = 300, height = 8, width = 16, units = ("in"), bg = 'white')
