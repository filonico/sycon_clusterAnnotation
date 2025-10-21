setwd("/scratch/evassvis/fn76/ANALYSIS/origin_of_vision/")

library(Seurat)
library(SeuratExtend)
library(reticulate)
library(scCustomize)
library(tidyverse)
library(metacell)


##############################
#     PLOT SAMAP RESULTS     #
##############################

# import the samap anndata object
ScilSlac_samap <- sceasy::convertFormat("05_SAMap/ScilSlac_leiden3Clusters_samap.h5ad", from = "anndata", to = "seurat")

# import the samap anndata object
AqueSlac_samap <- sceasy::convertFormat("05_SAMap/AqueSlac_leiden3Clusters_samap.h5ad", from = "anndata", to = "seurat")

# import the samap anndata object
AqueScil_samap <- sceasy::convertFormat("05_SAMap/AqueScil_leiden3Clusters_samap.h5ad", from = "anndata", to = "seurat")

# import Chris's gene pairs
ScilSlac_testPairs <- read.table("00_input/test_pairs.tsv", header = TRUE, sep = "\t") %>%
  # let IDs match the one from the anndata object
  mutate(Scil = paste0("Scil-", Scil),
         Slac = paste0("Slac-", Slac, "-g1"),
         # add an NA column to allow to plot just two features in FeaturePlot3.grid
         NA_column = NA,
         # add a column to match the SAMap gene pair file
         combined = paste0(Scil, "|", Slac)) %>%
  mutate(combined = str_remove(combined, "-g[0-9]$")) %>%
  relocate(NA_column, .after = 1)

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
