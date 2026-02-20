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

sycon_allMarkers_default %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5

heatmap <- DoHeatmap(Sycon, features = top5$gene) + NoLegend() +
  theme(axis.text = element_text(size = 5))
heatmap

top5 %>%
  write.table(file = "07_notableGenes_clusterAnnotation/top5_markers_perCluster.tsv",
              col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)


ggsave("07_notableGenes_clusterAnnotation/top5_marker_heatmap.pdf",
       heatmap, device = cairo_pdf,
       height = 8*2, width = 10*2, dpi = 300, bg = "white")
ggsave("07_notableGenes_clusterAnnotation/top5_marker_heatmap.png",
       heatmap, device = "png",
       height = 8*2, width = 10*2, dpi = 300, bg = "white")

cluster_21_26_specific <- sycon_allMarkers_default %>%
  filter(p_val_adj < 0.05) %>%
  group_by(gene) %>%
  mutate(cluster_agg = paste(unique(cluster), collapse = ",")) %>%
  ungroup() %>%
  filter(cluster_agg == "21,26") %>%
  slice_max(order_by = avg_log2FC, n = 50)

sycon_allMarkers_default %>%
  filter(p_val_adj < 0.05) %>%
  filter(gene == "g13150") %>%
  group_by(gene) %>%
  mutate(cluster_agg = paste(unique(cluster), collapse = ",")) %>%
  ungroup()
  
DoHeatmap(Sycon, features = cluster_21_26_specific$gene)  

uncertain_clusters <- cluster_identity %>% stack() %>%
  filter(str_detect(values, "Uncertain")) %>%
  pull(ind)

uncertain_clusters_specific #<-
  sycon_allMarkers_default %>%
  filter(p_val_adj < 0.05) %>%
  filter(cluster %in% c("1","7","11","16","21","23","25")) %>%
  group_by(gene) %>%
  mutate(n_clusters = n_distinct(cluster),
         cluster_agg = paste(unique(cluster), collapse = ",")) %>%
  ungroup() %>% View
  # filter(n_cl1+7+11+16+21+23+25usters > 5) %>% filter(gene == "g7730")
  # filter(cluster == "28",
  #        avg_log2FC > 1) #%>%
  View

DoHeatmap(Sycon, features = uncertain_clusters_specific$gene)  

cluster_5_specific <-
  sycon_allMarkers_default %>%
  filter(p_val_adj < 0.05) %>%
  group_by(gene) %>%
  mutate(n_clusters = n_distinct(cluster),
         cluster_agg = paste(unique(cluster), collapse = ",")) %>%
  ungroup() %>%
  filter(cluster == "5") %>%
  slice_max(order_by = avg_log2FC, n = 50)
  
DoHeatmap(Sycon, features = cluster_5_specific$gene)  

sycon_allMarkers_default %>%
  filter(p_val_adj < 0.05) %>%
  write.table(file = "07_notableGenes_clusterAnnotation/all_markers_perCluster.tsv",
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

# if a cluster significantly expresses known markers, we give it a name
# if a cluster maps with a spongilla cluster, we append the name 
# if a cluster does not express any known marker, we call it "unknown"
# if a cluster expresses known markers (but not significantly), we call it "uncertain"
# if a cluster significantly expresses known markers but with multiple identities, we call it "uncertain"
cluster_identity <- c("0" = "Uncertain_0", # uncertain identity (choanos + development)
                      "1" = "Choanocytes_1",
                      "2" = "Unknown_2",
                      "3" = "Uncertain_3", # not significant
                      "4" = "Uncertain_4", # not significant
                      "5" = "Uncertain_5", # uncertain identity (choanos + development)
                      "6" = "Unknown_6",
                      "7" = "Choanocytes_7",
                      "8" = "Uncertain_8", # not significant
                      "9" = "Uncertain_9", # not significant
                      "10" = "Uncertain_10", # uncertain identity (choanos + development)
                      "11" = "Choanocytes_11",
                      "12" = "Choanocytes_12",
                      "13" = "Choanocytes_13",
                      "14" = "Unknown_14",
                      "15" = "Metabolocytes-like_early_embryos",
                      "16" = "Choanocytes_16",
                      "17" = "Oocytes/early_embryos",
                      "18" = "Unknown_18",
                      "19" = "Oocytes/early_embryos",
                      "20" = "Unknown_20",
                      "21" = "Myopeptidocyte-like_choanocytes",
                      "22" = "Uncertain_22", # not significant
                      "23" = "Choanocytes_23",
                      "24" = "Oocytes/early_embryos",
                      "25" = "Sclerocytes-like_choanocytes",
                      "26" = "Uncertain_26", # uncertain identity (choanos + development)
                      "27" = "Oocytes/early_embryos",
                      "28" = "Uncertain_28", # not significant
                      "29" = "Choanocytes_29",
                      "30" = "Accessory_cells",
                      "31" = "Uncertain_31",
                      "32" = "Archaeocytes-like_stem_cells")
                      

as_tibble(cluster_identity) %>%
  rownames_to_column(var = "cluster_ID") %>%
  mutate(cluster_ID = as.character(as.numeric(cluster_ID) - 1)) %>%
  write.table("intermediate_results/01_sycon_Sobjects_umaps/cluster_identity.tsv", sep = "\t",
              quote = FALSE, col.names = TRUE, row.names = FALSE)

Sycon[[]] <- Sycon[[]] %>%
  left_join(data.frame(cluster_identity) %>%
              rownames_to_column(var = "cluster_n"),
            by = join_by("seurat_clusters" == "cluster_n"))

umap_clusters_raw <- DimPlot(Sycon, group.by = "seurat_clusters")
umap_clusters_raw

umap_clusters <- umap_clusters_raw@data %>%
  mutate(seurat_clusters = fct_infreq(as.factor(seurat_clusters))) %>%
  mutate(seurat_clusters = fct_rev(seurat_clusters)) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2, colour = seurat_clusters, label = seurat_clusters)) +
  geom_point(size = 1) +
  
  geom_label(data = . %>%
               group_by(seurat_clusters) %>%
               summarise(UMAP_1 = median(UMAP_1),
                         UMAP_2 = median(UMAP_2),
                         .groups = "drop"),
             aes(label = seurat_clusters),
             size = 3, linewidth = 0.5,
             text.color = "black", fill = alpha("white", 0.8)) +
  
  scale_color_manual(values = SeuratExtend::color_pro("default", n = 33,
                                                      sort = "diff")) + 
  
  labs(x = "UMAP 1", y = "UMAP 2",
       color = "Cell clusters") +
  
  annotation_custom(grob = segmentsGrob(x0 = unit(0, "mm"), x1 = unit(12, "mm"),
                                        y0 = unit(0, "mm"), y1 = unit(0, "mm"),
                                        arrow = arrow(length = unit(2.5, "mm"),
                                                      ends = "last", type = "open"),
                                        gp = gpar(col = "black", fill = "black", lwd = 1))) +
  
  annotation_custom(grob = segmentsGrob(x0 = unit(0, "mm"), x1 = unit(0, "mm"),
                                        y0 = unit(0, "mm"), y1 = unit(12, "mm"),
                                        arrow = arrow(length = unit(2.5, "mm"),
                                                      ends = "last", type = "open"),
                                        gp = gpar(col = "black", fill = "black", lwd = 1))) +
  
  annotation_custom(grob = textGrob(label = "UMAP 1",
                                    x = unit(0, "mm"), y = unit(0, "mm") - unit(2.5, "mm"),
                                    just = c(0, 1), gp = gpar(fontsize = 10))) +
  
  annotation_custom(grob = textGrob(label = "UMAP 2",
                                    x = unit(0, "mm") - unit(2.5, "mm"), y = unit(0, "mm"),
                                    just = c(0, 0), rot = 90, gp = gpar(fontsize = 10))) +
  
  coord_cartesian(clip = "off") +
  
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        plot.margin = margin(5, 5, 10, 8, "mm")) +
  theme_for_UMAPS
umap_clusters

umap_annotated_raw <- DimPlot(Sycon, group.by = "cluster_identity")

umap_annotated <- umap_annotated_raw@data %>%
  mutate(cluster_identity = str_replace_all(cluster_identity, "_", " "),
         cluster_identity = str_replace_all(cluster_identity, " [0-9]+$", ""),
         cluster_identity = str_replace(cluster_identity, "-like ", "-like\n"),
         cluster_identity = str_replace(cluster_identity, "/", "/\n")) %>%
  
  ggplot(aes(UMAP_1, UMAP_2)) +
  
  geom_point(data = . %>%
               subset(grepl("Unknown", cluster_identity)),
             colour = "grey85", size = 1) +
  geom_point(data = . %>%
               subset(grepl("Uncertain", cluster_identity)),
             colour = "grey40", size = 1) +
  geom_point(data = . %>%
               subset(!grepl("Unknown|Uncertain", cluster_identity)),
             aes(colour = as.factor(cluster_identity)), size = 1) +
  
  scale_color_manual(values = SeuratExtend::color_iwh("default", n = 7)) + 
  
  guides(color = guide_legend(override.aes = list(size = 2),
                              byrow = TRUE, position = "right",
                              keyheight = 1, default.unit = "cm")) +
  
  labs(colour = "Main cell types") +
  
  theme_bw(base_size = 12) +
  theme_for_UMAPS
umap_annotated

panel <- ggpubr::ggarrange(umap_clusters, umap_annotated, ncol = 2,
                           widths = c(0.8, 1))
panel

ggsave("07_notableGenes_clusterAnnotation/cluster_annotation_umaps.pdf",
       panel, device = cairo_pdf,
       width = 18/1.5, height = 8/1.5, units = "in", dpi = 300, bg = "white")
ggsave("07_notableGenes_clusterAnnotation/cluster_annotation_umaps.png",
       panel, device = "png",
       width = 18/1.5, height = 8/1.5, units = "in", dpi = 300, bg = "white")
