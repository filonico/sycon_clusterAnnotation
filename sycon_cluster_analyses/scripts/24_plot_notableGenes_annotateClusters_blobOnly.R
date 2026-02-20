#!/usr/bin/env Rscript

setwd("/lustre/alice3/data/evassvis/fn76/sycon/sycon_clusterAnnotation/sycon_cluster_analyses/")

library(Seurat)
library(DESeq2)
library(tidyverse)


#######################
#     LOAD INPUTS     #
#######################

Sycon_blobOnly <- readRDS("13_recluster_blob/Sycon_blobOnly.Rds")

sycon_blobOnly_allMarkers_default <- FindAllMarkers(Sycon_blobOnly,
                                                    group.by = "seurat_clusters_new",
                                                    only.pos = TRUE)

DimPlot(Sycon_blobOnly, group.by = "seurat_clusters_new")


###############################
#     DEFINE PLOT THEMES      #
###############################

theme_for_UMAPS <- theme(
  # aspect.ratio = 1,
  plot.background = element_blank(), 
  panel.border = element_blank(),
  panel.background = element_blank(),
  panel.grid = element_blank(),
  plot.title = element_text(size = 13, hjust = 0.5, vjust = 1.75, face = "bold"),
  axis.line = element_blank(),
  axis.ticks = element_blank(),
  axis.text = element_blank(),
  axis.title = element_blank()
)


#######################
#     GET MARKERS     #
#######################

sycon_blobOnly_allMarkers_default %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5

heatmap <- DoHeatmap(Sycon_blobOnly, features = top5$gene) + NoLegend() +
  theme(axis.text = element_text(size = 5))
heatmap

top5 %>%
  write.table(file = "13_recluster_blob/06_notableGenes_clusterAnnotation/top5_markers_perCluster.tsv",
              col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)


ggsave("13_recluster_blob/06_notableGenes_clusterAnnotation/top5_marker_heatmap.pdf",
       heatmap, device = cairo_pdf,
       height = 8*2, width = 10*2, dpi = 300, bg = "white")
ggsave("13_recluster_blob/06_notableGenes_clusterAnnotation/top5_marker_heatmap.png",
       heatmap, device = "png",
       height = 8*2, width = 10*2, dpi = 300, bg = "white")

sycon_blobOnly_allMarkers_default %>%
  filter(p_val_adj < 0.05) %>%
  write.table(file = "13_recluster_blob/06_notableGenes_clusterAnnotation/all_markers_perCluster.tsv",
              col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)


###############################
#     PLOT NOTABLE GENES      #
###############################

notable_genes <- list(
  "Sox" = c("g11067","g5646","g5684","g1695","g11976","g8388"),
  "Tgf-B\nligands" = c("g4396","g13016","g8882","g8564","g2075","g13405","g6992","g2931","g12087","g13365"),
  "Wnt" = c("g7201","g8624","g13150","g4399","g6515","g1564","g13140","g2924","g13308","g123","g6996","g12657","g13684","g7152","g11232","g1969"),
  "Pax/Six/\nEya" = c("g7040","g8356","g3655"),
  "Spicule\nformation" = c("g782","g11209","g10087","g8503","g7906","g9905", "g9130", "g6914", "g5755", "g5244", "g2825"),
  "Fzd\n(Wnt\npath)" = c("g4648","g3042"),
  "Dvl\n(Wnt\npath)" = c("g949"),
  "Lrp\n(Wnt\npath)" = c("g13476"),
  "Tcf\n(Wnt\npath)" = c("g9","g11378"),
  "Beta-cat\n(Wnt path)" = c("g13719","g4975"),
  "Smad" = c("g9893","g11243","g4903","g1511","g9154","g11177"),
  "Bra/Tbox" = c("g3032","g3118","g9948","g10519"),
  "Totipotency/\ngermline/\noocytes" = c("g3867","g9458", "g12083", "g11831", "g3684", "g12115", "g6724","g12152", "g3613"),
  "Elav" = c("g7460"),
  "Gata" = c("g2248"),
  "Gli" = c("g8905")
)

p <- Sycon_blobOnly %>% SeuratExtend::DotPlot2(features = notable_genes, flip = TRUE)

dotplot_notable_genes <- p$data %>%
  mutate(Var2 = factor(Var2,
                       levels = c("1", "7", "11", "12", "13", "16", "23", "21", "26", "29",
                                  "15", "25",
                                  "32", "27", "24", "19", "17",
                                  "22", "20", "28", "30", "31",
                                  "0", "2", "3", "4", "5", "6", "8", "9", "10", "14", "18"))) %>%
  filter(pct > 0) %>%
  # semi_join(sycon_allMarkers_default %>%
  #             filter(p_val_adj >= 0.05),
  #           by = c("Var1" = "gene", "Var2" = "cluster")) %>%
  semi_join(sycon_blobOnly_allMarkers_default,
            by = c("Var1" = "gene", "Var2" = "cluster")) %>%
  
  ggplot(aes(Var1, Var2)) +
  geom_point(aes(size = pct, fill = zscore), shape = 21) +
  
  scale_fill_viridis_c(option = "plasma") +
  
  scale_y_discrete(limits = rev) +
  scale_x_discrete(position = "top") +
  scale_size_continuous(range = c(2,9)) +
  
  xlab("Notable genes") +
  ylab("Cell clusters") +
  labs(size = "Percent\nexpressed") +
  
  facet_grid(cols = vars(FeatureGroup), scales = "free_x", space = "free") +
  
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 8),
        strip.clip = "off",
        axis.text.x = element_text(angle = 45, hjust = 0, size = 7),
        panel.border = element_rect(color = "#4f4f4f", fill = NA, linewidth = 0.8),
        panel.background = element_blank(),
        panel.grid = element_line(color = "#dbdbdb"))

dotplot_notable_genes

ggsave("13_recluster_blob/06_notableGenes_clusterAnnotation/dotplot_notableGenes.pdf",
       dotplot_notable_genes, device = cairo_pdf,
       dpi = 300, height = 8, width = 20, units = ("in"), bg = 'white')
ggsave("13_recluster_blob/06_notableGenes_clusterAnnotation/dotplot_notableGenes.png",
       dotplot_notable_genes, device = "png",
       dpi = 300, height = 8, width = 20, units = ("in"), bg = 'white')

dotplot_notable_genes_sig <- p$data %>%
  mutate(Var2 = factor(Var2,
                       levels = c("1", "7", "11", "12", "13", "16", "23", "21", "26", "29",
                                  "15", "25",
                                  "32", "27", "24", "19", "17",
                                  "22", "20", "28", "30", "31",
                                  "0", "2", "3", "4", "5", "6", "8", "9", "10", "14", "18"))) %>%
  filter(pct > 0) %>%
  semi_join(sycon_blobOnly_allMarkers_default %>%
              filter(p_val_adj >= 0.05),
            by = c("Var1" = "gene", "Var2" = "cluster")) %>%
  
  ggplot(aes(Var1, Var2)) +
  geom_point(aes(size = pct, fill = zscore), shape = 21) +
  
  scale_fill_viridis_c(option = "plasma") +
  
  scale_y_discrete(limits = rev) +
  scale_x_discrete(position = "top") +
  scale_size_continuous(range = c(2,9)) +
  
  xlab("Notable genes") +
  ylab("Cell clusters") +
  labs(size = "Percent\nexpressed") +
  
  facet_grid(cols = vars(FeatureGroup), scales = "free_x", space = "free") +
  
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 8),
        strip.clip = "off",
        axis.text.x = element_text(angle = 45, hjust = 0, size = 7),
        panel.border = element_rect(color = "#4f4f4f", fill = NA, linewidth = 0.8),
        panel.background = element_blank(),
        panel.grid = element_line(color = "#dbdbdb"))

dotplot_notable_genes_sig

ggsave("13_recluster_blob/06_notableGenes_clusterAnnotation/dotplot_notableGenes_sig.pdf",
       dotplot_notable_genes_sig, device = cairo_pdf,
       dpi = 300, height = 8, width = 15, units = ("in"), bg = 'white')
ggsave("13_recluster_blob/06_notableGenes_clusterAnnotation/dotplot_notableGenes_sig.png",
       dotplot_notable_genes_sig, device = "png",
       dpi = 300, height = 8, width = 15, units = ("in"), bg = 'white')