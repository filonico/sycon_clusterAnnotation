#!/usr/bin/env Rscript

setwd("/lustre/alice3/data/evassvis/fn76/sycon/sycon_clusterAnnotation/sycon_cluster_analyses/")

library(Seurat)
library(SeuratExtend)
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
  legend.text = element_text(size = 8),
  legend.title = element_text(size = 10, face = "bold"),
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

DotPlot2(Sycon, features = top5$gene)

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
  "Gli" = c("g8905"),
  "HD\ngenes" = c("g8932","g8162","g9066","g9084","g8884","g6644","g8884","g6644")
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
  semi_join(sycon_allMarkers_default,
            by = c("Var1" = "gene", "Var2" = "cluster")) %>%
  # left_join(notable_gene_names, by = join_by("Var1" == "gene_ID")) %>%
  
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
                      "3" = "Unknown_3",
                      "4" = "Unknown_4",
                      "5" = "Uncertain_5", # uncertain identity (choanos + development)
                      "6" = "Unknown_6",
                      "7" = "Choanocytes_7",
                      "8" = "Unknown_8",
                      "9" = "Unknown_9",
                      "10" = "Uncertain_10", # uncertain identity (choanos + development)
                      "11" = "Choanocytes_11",
                      "12" = "Unknown_12",
                      "13" = "Unknown_13",
                      "14" = "Unknown_14",
                      "15" = "Metabolocyte-_and_pinacocyte-like_early_embryos",
                      "16" = "Choanocytes_16",
                      "17" = "Oocytes/early_embryos",
                      "18" = "Unknown_18",
                      "19" = "Oocytes/early_embryos",
                      "20" = "Unknown_20",
                      "21" = "Myopeptidocyte-like_choanocytes",
                      "22" = "Unknown_22",
                      "23" = "Choanocytes_23",
                      "24" = "Oocytes/early_embryos",
                      "25" = "Sclerocyte-like_pinacocytes",
                      "26" = "Uncertain_26", # uncertain identity (choanos + development)
                      "27" = "Oocytes/early_embryos",
                      "28" = "Unknown_28",
                      "29" = "Unknown_29",
                      "30" = "Accessory_cells",
                      "31" = "Unknown_31",
                      "32" = "Archaeocyte-like_stem_cells")
                      

cluster_identity_df <- as_tibble(cluster_identity) %>%
  rownames_to_column(var = "cluster_ID") %>%
  mutate(cluster_ID = as.character(as.numeric(cluster_ID) - 1))

write.table(cluster_identity_df, "07_notableGenes_clusterAnnotation/cluster_identity.tsv", sep = "\t",
            quote = FALSE, col.names = TRUE, row.names = FALSE)

Sycon[[]] <- Sycon[[]] %>%
  left_join(cluster_identity_df,
            by = join_by("seurat_clusters" == "cluster_ID")) %>%
  rename("cluster_identity" = "value")

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
             size = 2.5, linewidth = 0.5,
             text.color = "black", fill = alpha("white", 0.8)) +
  
  scale_color_manual(values = SeuratExtend::color_pro("default", n = 33,
                                                      sort = "diff")) + 
  
  labs(x = "UMAP 1", y = "UMAP 2",
       color = "Cell clusters") +
  
  annotation_custom(grob = grid::segmentsGrob(x0 = unit(0, "mm"), x1 = unit(12, "mm"),
                                              y0 = unit(0, "mm"), y1 = unit(0, "mm"),
                                              arrow = arrow(length = unit(2.5, "mm"),
                                                            ends = "last", type = "open"),
                                              gp = grid::gpar(col = "black", fill = "black", lwd = 1))) +
  
  annotation_custom(grob = grid::segmentsGrob(x0 = unit(0, "mm"), x1 = unit(0, "mm"),
                                              y0 = unit(0, "mm"), y1 = unit(12, "mm"),
                                              arrow = arrow(length = unit(2.5, "mm"),
                                                            ends = "last", type = "open"),
                                              gp = grid::gpar(col = "black", fill = "black", lwd = 1))) +
  
  annotation_custom(grob = grid::textGrob(label = "UMAP 1",
                                          x = unit(0, "mm"), y = unit(0, "mm") - unit(2.5, "mm"),
                                          just = c(0, 1), gp = grid::gpar(fontsize = 10))) +
  
  annotation_custom(grob = grid::textGrob(label = "UMAP 2",
                                          x = unit(0, "mm") - unit(2.5, "mm"), y = unit(0, "mm"),
                                          just = c(0, 0), rot = 90, gp = grid::gpar(fontsize = 10))) +
  
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
         cluster_identity = str_replace(cluster_identity, "/", "/\n"),
         cluster_identity = str_replace(cluster_identity, "and ", "and \n")) %>%
  
  ggplot(aes(UMAP_1, UMAP_2)) +
  
  geom_point(data = . %>%
               subset(grepl("Unknown", cluster_identity)),
             aes(color = "Unknown"), size = 1) +
  geom_point(data = . %>%
               subset(grepl("Uncertain", cluster_identity)),
             aes(colour = "Uncertain"), size = 1) +
  geom_point(data = . %>%
               subset(!grepl("Unknown|Uncertain", cluster_identity)),
             aes(colour = as.factor(cluster_identity)), size = 1) +
  
  scale_color_manual(values = c(SeuratExtend::color_iwh("default", n = 7),
                                "grey60", "grey85")) + 
  
  guides(color = guide_legend(override.aes = list(size = 2),
                              byrow = TRUE, position = "right",
                              keyheight = 0.8, default.unit = "cm")) +
  
  labs(colour = "Main cell types") +
  
  theme_bw(base_size = 12) +
  theme_for_UMAPS
umap_annotated

panel_umaps <- ggpubr::ggarrange(umap_clusters, umap_annotated, ncol = 2,
                           widths = c(0.8, 1))
panel_umaps

ggsave("07_notableGenes_clusterAnnotation/cluster_annotation_umaps.pdf",
       panel_umaps, device = cairo_pdf,
       width = 18/1.5, height = 8/1.5, units = "in", dpi = 300, bg = "white")
ggsave("07_notableGenes_clusterAnnotation/cluster_annotation_umaps.png",
       panel_umaps, device = "png",
       width = 18/1.5, height = 8/1.5, units = "in", dpi = 300, bg = "white")


######################################
#     FIND MARKERS PER CELL TYPE     #
######################################

Sycon[[]] <- Sycon[[]] %>%
  mutate(cluster_identity_reduced = cluster_identity, 
         cluster_identity_reduced = str_replace_all(cluster_identity_reduced, "_", " "),
         cluster_identity_reduced = str_replace_all(cluster_identity_reduced, " [0-9]+$", ""))
  
sycon_allMarkers_cellTypes <- Sycon %>%
  FindAllMarkers(group.by = "cluster_identity_reduced")

sycon_allMarkers_cellTypes %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

heatmap <- DoHeatmap(Sycon,
                     group.by = "cluster_identity_reduced",
                     features = top10$gene) +
  NoLegend() +
  theme(axis.text = element_text(size = 5))
heatmap

top10 %>%
  filter(cluster == "Choanocytes") %>%
  pull("gene")


##############################################
#     PLOT NOTABLE GENES FOR PUBLICATION     #
##############################################

notable_gene_names <- data.frame(gene_name = c("AE-like1","βcatA","βcatB","Bra1","Bra2","CA1","CA2",
                                               "Cdx","Diactinin","DvlA","DvlB","Elav","Eya",
                                               "FzdA","FzdB","FzdC","FzdD","Gata","Gli","Hedgling",
                                               "Hex","Hmx","Lrp5/6","MsiA","MsiB","Msx",
                                               "Nanos","NCBT-like1","NKA","NKB","PaxB","Pl10A","Pl10B",
                                               "SixA","Smad1/5","Smad4","SmadRa",
                                               "Sox6","Sox7","SoxB","SoxC","SoxE","SoxF1","SoxF2","Spiculin",
                                               "TboxA","TboxB","TboxC","TcfA","TcfB",
                                               "TgfβB","TgfβD","TgfβF","TgfβG","TgfβN","Triactinin","VasaA","VasaB",
                                               "WntA","WntC","WntD","WntF","WntG","WntI","WntJ","WntL","WntN","WntQ","WntR","WntS","WntT"),
                                 gene_ID = c("g8503","g13719","g4975","g3032","g3118","g782","g11209","g8932",
                                             NA,"g949",NA,"g7460","g8356",NA,NA,"g4648","g3042","g2248","g8905",
                                             NA,"g8162","g9066","g13476","g11831","g12083","g9084","g3613","g10087",
                                             "g8884","g6644","g3655","g6724","g12152","g7040","g9893",
                                             "g11243","g9154","g11067","g5646","g5684","g1695","g11976","g8388",NA,
                                             "g7906","g9948",NA,"g10519","g9","g11378","g13016",NA,"g8564","g2075","g6992",
                                             "g9905","g3867","g9458","g7201","g13150","g4399",NA,NA,"g13140",NA,
                                             "g13308","g6996","g13684","g7152","g11232","g1969")) %>% as_tibble()


p <- Sycon %>% SeuratExtend::DotPlot2(features = notable_gene_names$gene_ID)

dotplot_markers_publication <- p$data %>%
  
  # sort y axis values to match cell types
  mutate(Var2 = factor(Var2, levels = c(
    "1","7","11","16","23","21","100",
    "15","101",
    "25","102",
    "17","19","24","27","32","103",
    "30","104",
    "0","3","4","5","8","9","10","12","13","22","26","28","29","31","105",
    "2","6","14","18","20"
  ))) %>%
  
  # remove genes with 0 expression in cells
  filter(pct > 0) %>%
  
  # add gene names
  left_join(notable_gene_names, by = join_by("Var1" == "gene_ID")) %>%
  filter(!is.na(gene_name)) %>%
  
  # remove genes that are not significant in any cell
  filter(! gene_name %in% c("SoxC", "SoxF1", "Lrp5/6", "WntD", "WntS", "MsiA", "Gata", "TcfA", "βcatA", "Bra2", "Pl10B")) %>%
  
  # add cell cluster names
  left_join(cluster_identity_df, by = join_by("group" == "cluster_ID")) %>%
  mutate(value = str_replace_all(value, "_", " "),
         value = str_replace_all(value, " [0-9]+$", ""),
         value = str_replace(value, "-like ", "-like\n"),
         value = str_replace(value, "/", "/\n"),
         value = str_replace(value, "and ", "and\n"),
         value = factor(value, levels = c("Choanocytes", "Oocytes/\nearly embryos", "Archaeocyte-like\nstem cells",
                                          "Metabolocyte- and\npinacocyte-like\nearly embryos",
                                          "Sclerocyte-like\npinacocytes", "Myopeptidocyte-like\nchoanocytes",
                                          "Accessory cells", "Uncertain", "Unknown"))) %>%
  rename("cell_identity" = "value") %>%
  
  # add p-values according to FindAllMarkers
  left_join(sycon_allMarkers_default %>%
              select(c(p_val_adj, cluster,gene)), by = join_by("group" == "cluster", "Var1" == "gene"),
            relationship = "many-to-many") %>%

  mutate(gene_name = factor(gene_name, levels = c("SoxE","Pl10A","Elav","WntQ","VasaA","VasaB","Smad1/5","Sox6","Sox7","SoxB","SoxC","SoxF1",
                                                  "TgfβB","TgfβF","TgfβG","TgfβN","WntA","WntC","WntD","WntI","WntL","WntN","WntR","WntS",
                                                  "WntT","Eya","SmadRa","MsiB","Gata","FzdC","FzdD","DvlA","TcfB","Cdx","Gli","PaxB",
                                                  "AE-like1","CA1","CA2","NCBT-like1","Spiculin","Triactinin","βcatB","TboxC","NKA","NKB",
                                                  "Nanos","Hex","Smad4","TboxA","Hmx","SixA","Bra1","MsiA","Lrp5/6","TcfA","βcatA","Bra2",
                                                  "Pl10B")),
         sig = case_when(p_val_adj < 0.05 ~ "*",
                         TRUE ~ NA)) %>%
  arrange(zscore) %>%
  
  
  ggplot(aes(gene_name, as_factor(Var2))) +
  
  geom_vline(xintercept = Inf, color = "#4f4f4f") +
  geom_rect(data = . %>% distinct(cell_identity),
            aes(fill = cell_identity), alpha = 0.2, inherit.aes = FALSE,
            xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  scale_fill_manual(values = c("Accessory cells" = "#B49440",
                               "Archaeocyte-like\nstem cells" = "#B65CBF",
                               "Choanocytes" = "#67A64E",
                               "Metabolocyte- and\npinacocyte-like\nearly embryos" = "#747CC9",
                               "Myopeptidocyte-like\nchoanocytes" = "#CA5F3E",
                               "Oocytes/\nearly embryos" = "#48B1A7",
                               "Sclerocyte-like\npinacocytes" = "#C8577B",
                               "Uncertain" = "grey60",
                               "Unknown" = NA),
                    guide = "none") +
  ggnewscale::new_scale_fill() +
  
  geom_point(aes(size = pct, color = zscore), shape = 19, alpha = 0.8) +
  geom_text(aes(label = sig), size = 3, color = "white", fontface = "bold", vjust = 0.65) +
  
  scale_color_distiller(palette = "Blues", direction = 1) +
  # scale_color_gradient(high = "#08306B", low = "#DEEBF7") +
  
  scale_y_discrete(limits = rev) +
  scale_x_discrete(position = "top") +
  scale_size_continuous(range = c(1,7)) +
  
  labs(x = "Marker genes", y = "Cell clusters",
       color = "z-score", size = "Percent expressed") +
  
  facet_wrap(. ~ cell_identity, space = "free_y", scales = "free_y", strip.position = "right") +
  
  coord_cartesian(clip = "off") +
  guides(size = guide_legend(label.position = "bottom", override.aes = list(color = "#77b3d6")),
         color = guide_colorbar(barheight = 1.4)) +
  
  theme_bw(base_size = 12) +
  theme(
    plot.background = element_rect(fill = "transparent", colour = NA), 
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", lineend = "round", linewidth = .5),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA),
    legend.key.width = unit(.4, "cm"),
    legend.key.height = unit(.4, "cm"),
    legend.position = "bottom",
    legend.title = element_text(size = 10, hjust = 0.5, face = "bold"),
    legend.title.position = "top",
    legend.text = element_text(size = 8),
    axis.line = element_blank(),
    axis.ticks = element_line(colour = "grey90", linewidth = .5),
    axis.ticks.length = unit(0.10, "cm"),
    axis.text.x = element_text(color = "black", angle = 45, face = "italic", hjust = 0,
                               margin = margin(t = 4, r = 0, b = 0, l = 0)),
    axis.text.y = element_text(color = "black",
                               margin = margin(t = 0, r = 4, b = 0, l = 0)),
    axis.title.y = element_text(angle = 90, size = 10, face = "bold",
                                margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.title.x = element_text(angle = 0, size = 10, face = "bold",
                                margin = margin(t = 10, r = 0, b = 0, l = 0)),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y = element_text(face = "bold", size = 8, angle = 0, hjust = 0),
    strip.clip = "off",
  )
dotplot_markers_publication

panel_umaps_noLegend <- ggpubr::ggarrange(umap_clusters,
                                          umap_annotated +
                                            ggrepel::geom_text_repel(data = data.frame(x = -2, y = -5.5, label = "Uncertain"),
                                                                     aes(x = x, y = y, label = label), nudge_x = -0.5, nudge_y = -2) +
                                            ggrepel::geom_text_repel(data = data.frame(x = 1, y = -3.4, label = "Unknown"),
                                                                     aes(x = x, y = y, label = label), nudge_x = 2.5, nudge_y = -2) +
                                            ggrepel::geom_text_repel(data = data.frame(x = 7.3, y = -10.4, label = "Myopeptidocyte-like\nchoanocytes"),
                                                                     aes(x = x, y = y, label = label), nudge_x = -1, nudge_y = 2.5) +
                                            ggrepel::geom_text_repel(data = data.frame(x = 9.4, y = -3.1, label = "Metabolocyte- and\npinacocyte-like\nearly-embryos"),
                                                                     aes(x = x, y = y, label = label), nudge_x = -4, nudge_y = 0.2, max.overlaps = Inf) +
                                            ggrepel::geom_text_repel(data = data.frame(x = 11.5, y = 0, label = "Sclerocyte-like\npinacocytes"),
                                                                     aes(x = x, y = y, label = label), nudge_x = -4, nudge_y = 0.1) +
                                            ggrepel::geom_text_repel(data = data.frame(x = 10.1, y = 1.9, label = "Accessory cells"),
                                                                     aes(x = x, y = y, label = label), nudge_x = -1.8, nudge_y = 1.2) +
                                            ggrepel::geom_text_repel(data = data.frame(x = 11.9, y = 5.2, label = "Archaeocyte-like\n stem cells"),
                                                                     aes(x = x, y = y, label = label), nudge_x = -4, nudge_y = 0.4) +
                                            ggrepel::geom_text_repel(data = data.frame(x = 7, y = 9.6, label = "Oocytes/early embryos"),
                                                                     aes(x = x, y = y, label = label), nudge_x = -2.5, nudge_y = -1) +
                                            ggrepel::geom_text_repel(data = data.frame(x = 0.9, y = 4.1, label = "Choanocytes"),
                                                                     aes(x = x, y = y, label = label), nudge_x = 2.5, nudge_y = 1.4) +
                                            NoLegend(), ncol = 2)
panel_umaps_noLegend

panel_fig1 <- ggpubr::ggarrange(panel_umaps_noLegend, dotplot_markers_publication,
                                nrow = 2, labels = "AUTO", heights = c(1,2))
panel_fig1

ggsave("07_notableGenes_clusterAnnotation/fig1_panel.pdf",
       panel_fig1, device = cairo_pdf,
       width = 20/1.5, height = 22/1.5, units = "in", dpi = 300, bg = "white")
ggsave("07_notableGenes_clusterAnnotation/fig1_panel.png",
       panel_fig1, device = "png",
       width = 20/1.5, height = 22/1.5, units = "in", dpi = 300, bg = "white")

