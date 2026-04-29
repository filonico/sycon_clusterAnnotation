#!/usr/bin/env Rscript

setwd("/data/evassvis/fn76/sycon/sycon_clusterAnnotation/sycon_cluster_analyses/")

library(tidyverse)
library(Seurat)
library(SeuratExtend)


######################
#     LOAD INPUT     #
######################

load("00_input/Sycon_Seuratv4.Rdata")

sycon_allMarkers_default <- FindAllMarkers(Sycon,
                                           only.pos = TRUE)

sycon_allMarkers_default %>%
  filter(cluster %in% c(1, 7, 11, 16, 23)) %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  select(gene) %>%
  distinct() %>%
  pull(gene)

choanocyte_markers <- c("g7939", "g13347", "g13351", "g13339", "g10860", "g11926", "g10859", "g4666", "g12872", "g3385", "g12871", "g6860", "g307", "g9627", "g8652")

SeuratExtend::DotPlot2(Sycon, features = choanocyte_markers)

SeuratExtend::DimPlot2(Sycon, features = choanocyte_markers)

###############################
#     DEFINE PLOT THEMES      #
###############################

umap_arrows <- list(
  
  annotation_custom(grob = grid::segmentsGrob(
    x0 = unit(0, "mm"), x1 = unit(12, "mm"),
    y0 = unit(0, "mm"), y1 = unit(0, "mm"),
    arrow = arrow(length = unit(2.5, "mm"), ends = "last", type = "open"),
    gp = grid::gpar(col = "black", fill = "black", lwd = 1))),
  
  annotation_custom(grob = grid::segmentsGrob(
    x0 = unit(0, "mm"), x1 = unit(0, "mm"),
    y0 = unit(0, "mm"), y1 = unit(12, "mm"),
    arrow = arrow(length = unit(2.5, "mm"), ends = "last", type = "open"),
    gp = grid::gpar(col = "black", fill = "black", lwd = 1))),
  
  annotation_custom(grob = grid::textGrob(label = "UMAP 1",
                                          x = unit(0, "mm"), y = unit(0, "mm") - unit(2.5, "mm"),
                                          just = c(0, 1), gp = grid::gpar(fontsize = 10))),
  
  annotation_custom(grob = grid::textGrob(label = "UMAP 2",
                                          x = unit(0, "mm") - unit(2.5, "mm"), y = unit(0, "mm"),
                                          just = c(0, 0), rot = 90, gp = grid::gpar(fontsize = 10))),
  
  coord_cartesian(clip = "off"))

umap_coord_x <- list(annotation_custom(grob = grid::textGrob(label = "UMAP 1",
                                                             x = unit(0, "mm"), y = unit(0, "mm") - unit(2.5, "mm"),
                                                             just = c(0, 1), gp = grid::gpar(fontsize = 10))),
                     coord_cartesian(clip = "off"))

umap_coord_y <- list(annotation_custom(grob = grid::textGrob(label = "UMAP 2",
                                                             x = unit(0, "mm") - unit(2.5, "mm"), y = unit(0, "mm"),
                                                             just = c(0, 0), rot = 90, gp = grid::gpar(fontsize = 10))),
                     coord_cartesian(clip = "off"))

theme_for_UMAPS <- theme(
  # aspect.ratio = 1,
  plot.background = element_blank(), 
  panel.border = element_rect(color = "#545454"),
  panel.grid = element_blank(),
  legend.position = "none",
  plot.title = element_text(size = 13, hjust = 0, vjust = 1.75, face = "bold.italic"),
  axis.line = element_blank(),
  axis.ticks = element_blank(),
  axis.text = element_text(color = "#545454"),
  axis.title = element_blank()
)


#####################
#     FUNCTIONS     #
#####################

plot_dimplot_gene <- function(gene_id, gene_name) {
  
  dimplot_raw <- Sycon %>%
    Seurat::FeaturePlot(features = gene_id)
  
  dimplot_def <- dimplot_raw@data %>%
    arrange(.data[[gene_id]]) %>%
    ggplot(aes(UMAP_1, UMAP_2)) +
    geom_vline(xintercept = seq(-5, 10, 5), col = "grey95") +
    geom_hline(yintercept = seq(-10, 10, 5), col = "grey95") +
    geom_point(aes(col = .data[[gene_id]]), size = 1) +
    
    # scale_color_distiller(palette = "Blues", direction = 1) +
    scale_color_gradient(high = "darkblue", low = "grey80") +
    
    scale_y_continuous(breaks = seq(-10, 5, 5)) +
    
    annotate("text", x = -7.3, y = 9.9, hjust = 0,
             label = paste0("italic('", gene_name, "')~'(", gene_id, ")'"),
             parse = TRUE) +
    
    # labs(title = paste0(gene_name, " (", gene_id, ")")) +
    
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme_for_UMAPS
  
  return(dimplot_def)
}


############################
#     PLOT HCR MARKERS     #
############################

hcr_markers <- c(
  # "Vasa-B" = "g9458",
  # "Elav" = "g7460",
  "Granulin" = "g11545",
  # "Calcineurin" = "g2493",
  # "Actin" = "g13004",
  "Wnt-N" = "g6996",
  "TgfB-F" = "g8564",
  # "Reeler" = "g5618",
  # "ApoL" = "g7154",
  "Mpeg-1" = "g4312",
  "Spiculin" = "g7906",
  "Triactinin" = "g9905",
  # "Tbox-A" = "g9948",
  "Laminin" = "g7679",
  "Sox-E" = "g11976"
  )


dimplots_hcr_markers <- Map(plot_dimplot_gene, unname(hcr_markers), names(hcr_markers))

dimplots_hcr_markers$g7906 <- dimplots_hcr_markers$g7906 +
  scale_x_continuous(breaks = seq(0, 10, 5)) +
  scale_y_continuous(breaks = seq(-5, 5, 5)) +
  umap_coord_x + umap_coord_y

panel <- patchwork::wrap_plots(dimplots_hcr_markers, nrow = 2) +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 14, face = "bold"),
        plot.tag.location = "plot")
panel

ggsave("07_notableGenes_clusterAnnotation/fig2_hcrMarkers.png",
       panel, device = "png",
       height = 8/1.5, width = 17/1.5, units = "in", dpi = 300, bg = "white")
ggsave("07_notableGenes_clusterAnnotation/fig2_hcrMarkers.pdf",
       panel, device = cairo_pdf,
       height = 8/1.5, width = 17/1.5, units = "in", dpi = 300, bg = "white")
