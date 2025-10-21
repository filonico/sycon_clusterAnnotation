library(Seurat)
library(SeuratExtend)
library(tidyverse)
library(ggvenn)

setwd("/lustre/alice3/data/evassvis/fn76/sycon/sycon_clusterAnnotation/")

#####################
#     FUNCTIONS     #
#####################

get_markers_and_save <- function(seurat_object, group, output_dir) {
  markers <- seurat_object %>% FindAllMarkers(group.by = group, only.pos = TRUE)
  
  markers_list <- markers %>%
    filter(avg_log2FC > 0) %>%
    group_by(cluster) %>%
    # summarise(count = n())
    summarise(genes = list(gene), .groups = "drop") %>%
    { setNames(.$genes, .$cluster) }
  
  # write the list of upregulated genes per cluster to a file
  for (cluster_name in names(markers_list)) {
    file_path <- file.path(output_dir,
                           paste0("cluster", cluster_name, "_upregulatedGenes.ls"))
    writeLines(markers_list[[cluster_name]],
               file_path)
  }
  
  return(markers)
}

get_gene_universe <- function(seurat_object, out_filename){
  
  seurat_object@assays$RNA$counts %>%
    as_tibble(rownames = NA) %>%
    rownames_to_column(var = "gene") %>%
    mutate(sum = rowSums(across(where(is.numeric)))) %>%
    filter(sum > 0) %>%
    select(c("gene")) %>%
    
    write.table(file = out_filename,
                col.names = FALSE, row.names = FALSE, quote = FALSE)
}


#################
#     INPUT     #
#################

load("00_input/Sycon_Seuratv4.Rdata")


######################
#     RE-CLUSTER     #
######################

DefaultAssay(Sycon) <- "RNA"

Sycon %>%
  SeuratExtend::DimPlot2(group.by = "seurat_clusters", label = TRUE, repel = TRUE,
                         theme = theme_umap_arrows())

clusters_to_remove <- c("15", "17", "19", "21", "24", "26", "27", "25", "30", "32")

Sycon_blobOnly <- Sycon %>%
  UpdateSeuratObject() %>%
  DietSeurat(assays = "RNA") %>%
  subset(!(seurat_clusters %in% clusters_to_remove)) %>%
  SCTransform()

Sycon_blobOnly <- Sycon_blobOnly %>%
  RunPCA()

ElbowPlot(Sycon_blobOnly, ndims = 50) %>%
  ggsave(plot = ., filename = "13_recluster_blob/elbowplot.pdf", device = cairo_pdf)

Sycon_blobOnly <- Sycon_blobOnly %>%
  RunUMAP(dims = 1:30)

Sycon_blobOnly <- Sycon_blobOnly %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters()

Sycon_blobOnly %>% DimPlot()

Sycon_blobOnly@meta.data <- Sycon_blobOnly@meta.data %>%
  rownames_to_column() %>%
  left_join(Sycon@meta.data %>%
              rownames_to_column() %>%
              select(c("rowname", "seurat_clusters")),
            by = c("rowname" = "rowname"),
            # keep = FALSE,
            suffix = c("_new", "_original")) %>%
  column_to_rownames(var = "rowname")

saveRDS(Sycon_blobOnly, file = "13_recluster_blob/Sycon_blobOnly.Rds")
Sycon_blobOnly <- readRDS("13_recluster_blob/Sycon_blobOnly.Rds")

scCustomize::as.anndata(x = Sycon_blobOnly, main_layer = "counts",
                        other_layers = NULL, file_path = "13_recluster_blob/04_SAMap/",
                        file_name = "Scil_blobOnly.h5ad")


#####################
#     DIM PLOTS     #
#####################

# dimplot with the original cluster names
dimplot_clusters_original <- Sycon_blobOnly %>%
  SeuratExtend::DimPlot2(group.by = "seurat_clusters_original", pt.size = 1, cols = "auto",
                         label = TRUE, repel = TRUE, box = TRUE, label.color = "black",
                         theme = list(labs(title = "Original clusters"),
                                      theme_classic(),
                                      theme_umap_arrows()))
  
dimplot_clusters_original
  
# dimplot with the re-clustered cluster names
dimplot_clusters_new <- Sycon_blobOnly %>%
  SeuratExtend::DimPlot2(group.by = "seurat_clusters_new", pt.size = 1,
                         label = TRUE, repel = TRUE, box = TRUE, label.color = "black",
                         theme = list(labs(title = "New clusters"),
                                      theme_classic(),
                                      NoAxes()))
dimplot_clusters_new

# get the panel
panel_clusters <- ggpubr::ggarrange(dimplot_clusters_original, dimplot_clusters_new,
                  ncol = 2)

panel_clusters

ggsave(plot = panel_clusters, filename = "13_recluster_blob/dimplot_blobOnly_clusters_panel.pdf",
       device = cairo_pdf, dpi = 300, height = 6, width = 12, units = ("in"), bg = 'white')
ggsave(plot = panel_clusters, filename = "13_recluster_blob/dimplot_blobOnly_clusters_panel.png",
       device = "png", dpi = 300, height = 6, width = 12, units = ("in"), bg = 'white')

# get the new cluster composition relative to the original ones
new_cluster_composition <- Sycon_blobOnly[[]] %>%
  select(seurat_clusters_new, seurat_clusters_original) %>%
  group_by(seurat_clusters_new, seurat_clusters_original) %>%
  count()

# plot the cluster composition
cluster_composition_seurat_style <- ClusterDistrBar(Sycon_blobOnly$seurat_clusters_new, Sycon_blobOnly$seurat_clusters_original, percent = FALSE)
cluster_composition_seurat_style

color_palette <- ggplot_build(cluster_composition_seurat_style)$data[[1]]$fill %>% unique()

new_cluster_composition_col <- new_cluster_composition %>%
  ggplot(aes(x = n, y = seurat_clusters_new, fill = seurat_clusters_original)) +
  geom_col(position = "stack", colour = "white") +
  
  labs(#title = "New cluster composition",
       fill = "Original clusters") +
  ylab("New clusters") +
  xlab("Number of cells") +
  
  scale_fill_manual(values = c(color_palette)) +
  scale_x_continuous(expand = c(0, 0)) +
  
  theme_minimal() +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line = element_line(linewidth = 0.5),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 9, face = "bold"),
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 8))

new_cluster_composition_col

# ggsave("13_recluster_blob/new_cluster_composition_barplot.pdf",
#        new_cluster_composition_col, device = cairo_pdf,
#        dpi = 300, height = 6, width = 12, units = ("in"), bg = 'white')

new_cluster_composition_tile <- table(Sycon_blobOnly$seurat_clusters_original, Sycon_blobOnly$seurat_clusters_new) %>%
  as.data.frame() %>%
  rename(original_clusters = Var1, new_clusters = Var2, count = Freq) %>%
  group_by(new_clusters) %>%
  mutate(normalized = count / sum(count)) %>%  # proportion of cluster
  ungroup() %>%
  filter(count != 0) %>%
  
  ggplot(aes(new_clusters, original_clusters, fill = normalized)) +
  geom_tile() +
  geom_text(aes(label = count), color = "black", size = 3) +
  
  scale_fill_gradient(low = alpha("white", 0), high = "#37d461", limits = c(0,1),
                      breaks = seq(0, 1, 0.25), labels = format(seq(0, 1, 0.25), nsmall = 2)) +
  labs(#title = "New cluster composition",
       fill = "Normalized\nvalues") +
  xlab("New clusters") +
  ylab("Original clusters") +
  coord_flip() +
  
  theme_minimal() +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line = element_line(linewidth = 0.5),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 9, face = "bold"),
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 8))

new_cluster_composition_tile

# ggsave("13_recluster_blob/new_cluster_composition_tile.pdf",
#        new_cluster_composition_tile, device = cairo_pdf,
#        dpi = 300, height = 7, width = 9, units = ("in"), bg = 'white')

cluster_composition_panel <- ggpubr::ggarrange(new_cluster_composition_col, new_cluster_composition_tile,
                  ncol = 2, align = "hv", widths = c(2,1.5)) %>%
  ggpubr::annotate_figure(top = ggpubr::text_grob("Cluster composition", face = "bold"))

ggsave("13_recluster_blob/new_cluster_composition_panel.pdf",
       cluster_composition_panel, device = cairo_pdf,
       dpi = 300, height = 7, width = 18, units = ("in"), bg = 'white')
ggsave("13_recluster_blob/new_cluster_composition_panel.png",
       cluster_composition_panel, device = "png",
       dpi = 300, height = 7, width = 18, units = ("in"), bg = 'white')


##################################
#     GET MARKERS PER CLUSTER    #
##################################

markers_blobOnly_originalClusters <- get_markers_and_save(Sycon_blobOnly,
                                                          "seurat_clusters_original",
                                                          "13_recluster_blob/01_onlyBlob_originalClusters/")

markers_blobOnly_newClusters <- get_markers_and_save(Sycon_blobOnly,
                                                     "seurat_clusters_new",
                                                     "13_recluster_blob/02_onlyBlob_newClusters/")

markers_venn <- ggvenn::ggvenn(list("Original cluster\nmarkers" = markers_blobOnly_originalClusters$gene,
                                    "New cluster\nmarkers" = markers_blobOnly_newClusters$gene),
                               auto_scale = TRUE)

ggsave("13_recluster_blob/venn_markers.pdf",
       markers_venn, device = cairo_pdf,
       dpi = 300, height = 6, width = 6, units = ("in"), bg = 'white')
ggsave("13_recluster_blob/venn_markers.png",
       markers_venn, device = "png",
       dpi = 300, height = 6, width = 6, units = ("in"), bg = 'white')

#############################
#     GET GENE UNIVERSE     #
#############################

get_gene_universe(Sycon_blobOnly, "13_recluster_blob/sycon_onlyBlob_geneUniverse.ls")
