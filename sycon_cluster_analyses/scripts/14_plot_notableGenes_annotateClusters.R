#!/usr/bin/env Rscript

setwd("/lustre/alice3/data/evassvis/fn76/sycon/sycon_clusterAnnotation/sycon_cluster_analyses/")

library(Seurat)
library(DESeq2)
library(tidyverse)


#######################
#     LOAD INPUTS     #
#######################

load("00_input/Sycon_Seuratv4.Rdata")

sycon_allMarkers_default <- FindAllMarkers(Sycon,
                                           only.pos = TRUE)


#######################
#     GET MARKERS     #
#######################

# sycon_allMarkers_roc <- FindAllMarkers(Sycon,
#                                        test.use = "roc",
#                                        only.pos = TRUE)
# 
# sycon_allMarkers_roc %>%
#   group_by(cluster) %>%
#   dplyr::filter(avg_log2FC > 1) %>%
#   slice_head(n = 5) %>%
#   ungroup() -> top5_roc

sycon_allMarkers_default %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5

heatmap <- DoHeatmap(Sycon, features = top5$gene) + NoLegend() +
  theme(axis.text = element_text(size = 5))
heatmap

ggsave("07_notableGenes_clusterAnnotation/top5_marker_heatmap.pdf",
       heatmap, device = cairo_pdf,
       height = 8*2, width = 10*2, dpi = 300, bg = "white")
ggsave("07_notableGenes_clusterAnnotation/top5_marker_heatmap.png",
       heatmap, device = "png",
       height = 8*2, width = 10*2, dpi = 300, bg = "white")

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

p <- Sycon %>% SeuratExtend::DotPlot2(features = notable_genes, flip = TRUE)

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
  semi_join(sycon_allMarkers_default,
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

ggsave("07_notableGenes_clusterAnnotation/dotplot_notableGenes.pdf",
       dotplot_notable_genes, device = cairo_pdf,
       dpi = 300, height = 8, width = 20, units = ("in"), bg = 'white')
ggsave("07_notableGenes_clusterAnnotation/dotplot_notableGenes.png",
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
  semi_join(sycon_allMarkers_default %>%
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

ggsave("07_notableGenes_clusterAnnotation/dotplot_notableGenes_sig.pdf",
       dotplot_notable_genes_sig, device = cairo_pdf,
       dpi = 300, height = 8, width = 15, units = ("in"), bg = 'white')
ggsave("07_notableGenes_clusterAnnotation/dotplot_notableGenes_sig.png",
       dotplot_notable_genes_sig, device = "png",
       dpi = 300, height = 8, width = 15, units = ("in"), bg = 'white')


##################################
#     PLOT ANNOTATED_MARKERS     #
##################################


cluster_identity <- c("0" = "Uncertain_0",
                      "1" = "Unknown_1",
                      "2" = "Uncertain_2",
                      "3" = "Uncertain_3",
                      "4" = "Uncertain_4",
                      "5" = "Choanocytes-like_1",
                      "6" = "Unknown_6",
                      "7" = "Uncertain_7",
                      "8" = "Uncertain_8",
                      "9" = "Uncertain_9",
                      "10" = "Uncertain_10",
                      "11" = "Unknown_11",
                      "12" = "Choanocytes-like_2",
                      "13" = "Choanocytes-like_2",
                      "14" = "Unknown_14",
                      "15" = "Metabolocyte-like_choanocytes",
                      "16" = "Unknown_16",
                      "17" = "Oocytes/early_embryos",
                      "18" = "Unknown_18",
                      "19" = "Oocytes/early_embryos",
                      "20" = "Uncertain_20",
                      "21" = "Myopeptidocyte-like_choanocytes",
                      "22" = "Uncertain_22",
                      "23" = "Unknown_23",
                      "24" = "Oocytes/early_embryos",
                      "25" = "Sclerocytes",
                      "26" = "Choanocytes-like_2",
                      "27" = "Oocytes/early_embryos",
                      "28" = "Uncertain_28",
                      "29" = "Choanocytes-like_2",
                      "30" = "Choanocytes-like_pinacocytes-like",
                      "31" = "Uncertain_31",
                      "32" = "Oocytes/early_embryos")

as_tibble(cluster_identity) %>%
  rownames_to_column(var = "cluster_ID") %>%
  mutate(cluster_ID = as.character(as.numeric(cluster_ID) - 1)) %>%
  write.table("intermediate_results/01_sycon_Sobjects_umaps/cluster_identity.tsv", sep = "\t",
              quote = FALSE, col.names = TRUE, row.names = FALSE)

Sycon[[]] <- Sycon[[]] %>%
  left_join(data.frame(cluster_identity) %>%
              rownames_to_column(var = "cluster_n"),
            by = join_by("seurat_clusters" == "cluster_n"))

umap <- Sycon %>%
  SetIdent(value = "seurat_clusters") %>%
  DimPlot2(pt.size = 2, cols = "light",
           label = TRUE, label.size = 3, label.color = "black", box = TRUE,
           theme = list(labs(title = expression(paste(bolditalic("Sycon ciliatium"), bold(" cell atlas")))),
                        theme_classic(), theme_umap_arrows(), NoLegend()))
umap

umap_annotated <- DimPlot(Sycon, group.by = "cluster_identity") +
  theme_umap_arrows()
umap_annotated

umap_cell_types <- umap_annotated@data %>%
  ggplot(aes(UMAP_1, UMAP_2)) +
  
  geom_point(data = subset(umap_annotated@data, grepl("Unknown", cluster_identity)),
             colour = "grey85", size = 2) +
  geom_point(data = subset(umap_annotated@data, grepl("Uncertain", cluster_identity)),
             colour = "grey40", size = 2) +
  geom_point(data = subset(umap_annotated@data, !grepl("Unknown|Uncertain", cluster_identity)),
             aes(colour = as.factor(cluster_identity)), size = 2) +
  
  # scale_colour_discrete(name = "",
  labs(colour = "Main cell types") +
  
  theme_classic() +
  theme_umap_arrows()
umap_cell_types

panel <- ggpubr::ggarrange(umap, umap_cell_types, ncol = 2, align= "hv")

ggsave("07_notableGenes_clusterAnnotation/cluster_annotation_umaps.pdf",
       panel, device = cairo_pdf,
       width = 20, height = 8, units = "in", dpi = 300, bg = "white")
ggsave("07_notableGenes_clusterAnnotation/cluster_annotation_umaps.png",
       panel, device = "png",
       width = 20, height = 8, units = "in", dpi = 300, bg = "white")
