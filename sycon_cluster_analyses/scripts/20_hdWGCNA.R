#!/usr/bin/env Rscript

setwd("/data/evassvis/fn76/sycon/sycon_clusterAnnotation/sycon_cluster_analyses/")

suppressMessages({library(tidyverse)
  library(hdWGCNA)
  library(SeuratExtend)})


######################
#     LOAD INPUT     #
######################

# load("00_input/Sycon_Seuratv4.Rdata")
# 
# Sycon <- UpdateSeuratObject(Sycon)

Sycon <- readRDS("12_hdWGCNA/Sycon_hdWGCNA.Rds")


#####################
#     FUNCTIONS     #
#####################

plot_module_UMAP <- function(module_name) {
  
  # first get the coordinates of points in the umap
  umap <- Sycon@reductions$umap@cell.embeddings %>%
    as_tibble(rownames = NA) %>%
    rownames_to_column() %>%
    
    # then add information from meta data, including the score for each module
    left_join(Sycon@meta.data %>% rownames_to_column()) %>%
    select(rowname, UMAP_1, UMAP_2, all_of(module_name)) %>%
    
    # pivot data longer
    pivot_longer(-c(rowname, UMAP_1, UMAP_2), names_to = "module", values_to = "expr") %>%
    arrange(expr) %>%
    # group_by(module) %>%
    
    # plot data
    ggplot(aes(UMAP_1, UMAP_2)) +
    geom_point(aes(col = expr),
               # alpha = 0.5, 
               size = 0.7
    ) +
    
    # define color gradient rules
    # note that in the original hdWGCNA script they make the range symmetrical
    # is it legit? Isn't this affecting visualization?
    # here, we do not make the range symmetrical
    scale_color_gradient2(low = "grey75", mid = "grey95", high = module_name,
                          midpoint = 0) +
    
    ggtitle(module_name) +
    
    # make the plot squared
    coord_fixed() +
    
    theme(axis.text = element_text(size = 7),
          axis.title = element_text(size = 8),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          panel.border = element_rect(color = "#4f4f4f", fill = NA, linewidth = 0.8),
          panel.background = element_blank(),
          panel.grid = element_blank()
    )
  
  return(umap)
}


##########################################
#     hdWGCNA with standard pipeline     #
##########################################

DefaultAssay(Sycon) <- "SCT"
Sycon

# standard hdWGCNA pipeline
# pipeline with variable features is better than the other ("fraction")
# it retrieves better defined clusters
Sycon <- SetupForWGCNA(Sycon,
                       gene_select = "variable",
                       wgcna_name = "RNA_variableFeatures")

# Sycon <- SetupForWGCNA(Sycon,
#                        gene_select = "fraction",
#                        fraction = 0.05,
#                        wgcna_name = "RNA_fraction")

Sycon@misc$active_wgcna <- "RNA_variableFeatures"

# compute metacells by splitting both by samples and seurat clusters
# reduce number of cells to consider a cluster and the number of overlapping cells to 10
# otherwise it get rids of cluster 28..32
Sycon <- MetacellsByGroups(seurat_obj = Sycon,
                           group.by = c("orig.ident", "seurat_clusters"),
                           ident.group = "seurat_clusters",
                           min_cells = 10, k = 10,
                           assay = "RNA",
                           wgcna_name = "RNA_variableFeatures")

Sycon <- NormalizeMetacells(Sycon)

# it is important to use the RNA assay for DGE analysis, in order to use raw counts
# instead of normalized and/or integrated matrices
Sycon <- SetDatExpr(Sycon,
                    group_name = levels(Sycon[[]]$seurat_clusters),
                    group.by = "seurat_clusters",
                    assay = "RNA",
                    layer = "data")

# test different soft powers
Sycon <- TestSoftPowers(Sycon)

# # plot the results
# plotSoftPower_list <- PlotSoftPowers(Sycon)
# 
# # assemble with patchwork
# patchwork::wrap_plots(plotSoftPower_list, ncol = 2)

# construct co-expression network
Sycon <- ConstructNetwork(Sycon,
                          tom_outdir = "12_hdWGCNA/01_RNA_assay/",
                          tom_name = "RNA_variableFeatures", overwrite_tom = TRUE)

# PlotDendrogram(Sycon)

# this is necessary to include gorup by samples  in the next step
Sycon <- ScaleData(Sycon,
                   features = VariableFeatures(Sycon))

# compute all MEs in the full single-cell dataset
Sycon <- ModuleEigengenes(Sycon,
                          group.by.vars = "orig.ident")

# compute eigengene-based connectivity (kME):
Sycon <- ModuleConnectivity(Sycon,
                            group.by = "seurat_clusters")

# make a featureplot of hMEs for each module
# here, plotted ranges are symmetrical (as hdWGCNA default settings)
RNA_featurePlot_list <- ModuleFeaturePlot(Sycon,
                                          features = "hMEs",
                                          order = TRUE)

# plot patchworked plots
patchwork::wrap_plots(RNA_featurePlot_list,
                      ncol = 3)

# extract hME values and add them to meta data
hMEs <- Sycon %>%
  GetMEs(harmonized = TRUE)
Sycon@meta.data <- cbind(Sycon@meta.data, hMEs)

# create a custom UMAP for each module
# here, plotted ranges are not corrected, thus not symmetrical
RNA_umap_list <- list()
for (i in names(RNA_featurePlot_list)) {
  RNA_umap_list[[i]] <- plot_module_UMAP(i)
}

# plot patchworked plots
umaps <- print(patchwork::wrap_plots(RNA_umap_list, ncol = 3))
umaps

ggsave("12_hdWGCNA/01_RNA_assay/gene_module_umaps.pdf",
       umaps,  device = cairo_pdf,
       dpi = 300, height = 6, width = 12, units = ("in"), bg = 'white')

ggsave("12_hdWGCNA/01_RNA_assay/gene_module_umaps.png",
       umaps,  device = "png",
       dpi = 300, height = 6, width = 12, units = ("in"), bg = 'white')

# get the list of genes per module
gene_list <- Sycon@misc$RNA_variableFeatures$wgcna_modules %>%
  group_by(module) %>%
  summarise(genes = list(gene_name), .groups = "drop") %>%
  { setNames(.$genes, .$module) }

# write the list of genes per module to separate files
for (module in names(gene_list)) {
  file_path <- file.path("12_hdWGCNA/01_RNA_assay/", paste0(module, "_RNA_module_genes.ls"))
  writeLines(gene_list[[module]],
             file_path)
}

saveRDS(Sycon, "12_hdWGCNA/Sycon_hdWGCNA.Rds")

RNA_modules <- GetModules(Sycon)
RNA_mods <- levels(RNA_modules$module); mods <- mods[mods != 'grey']

RNA_dotplot <- Sycon %>%
  DotPlot(features = RNA_mods, group.by = 'seurat_clusters')


RNA_dotplot_custom <- RNA_dotplot$data %>%
  mutate(id = factor(id,
                     levels = c("32", "27", "24", "19", "17",
                                "30", "25", "15",
                                "26", "21",
                                "23", "18", "16", "14", "11",
                                "31", "29", "28", "22", "20", "13", "12",
                                seq(10, 0, -1))),
         features.plot = factor(features.plot,
                                levels = c("yellow", "turquoise", "brown", "green",
                                           "blue", "red"))) %>%
  filter(pct.exp > 0) %>%
  ggplot(aes(features.plot, id)) +
  geom_point(aes(size = pct.exp, fill = avg.exp.scaled), shape = 21) +
  
  scale_y_discrete(limits = rev) +
  scale_size_continuous(range = c(1,7)) +
  
  xlab("gene co-expression\nmodules") +
  ylab("Cell clusters") +
  labs(size = "Percent\nexpressed", fill = "Average\nexpression") +
  
  scale_fill_viridis_c(option = "plasma") +
  
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 8),
        strip.clip = "off",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        panel.border = element_rect(color = "#4f4f4f", fill = NA, linewidth = 0.8),
        panel.background = element_blank(),
        panel.grid = element_line(color = "#dbdbdb"))

RNA_dotplot_custom

ggsave("12_hdWGCNA/01_RNA_assay/RNA_gene_module_dotplot.pdf",
       RNA_dotplot_custom,  device = cairo_pdf,
       dpi = 300, height = 6, width = 12, units = ("in"), bg = 'white')

ggsave("12_hdWGCNA/01_RNA_assay/RNA_gene_module_dotplot.png",
       RNA_dotplot_custom,  device = "png",
       dpi = 300, height = 6, width = 12, units = ("in"), bg = 'white')


# ##########################################
# #     hdWGCNA with SCTransformed data    #
# ##########################################
# 
# # THE DEVELOPERS DO NOT SUGGEST USING SCTransformed DATA
# DefaultAssay(Sycon) <- "SCT"
# 
# Sycon <- Sycon %>%
#   SetupForWGCNA(features = VariableFeatures(Sycon),
#                 wgcna_name = "SCT_variableFeature")
# 
# Sycon@misc$active_wgcna <- "SCT_variableFeature"
# 
# # compute metacells by splitting both by samples and seurat clusters
# # reduce number of cells to consider a cluster and the number of overlapping cells to 10
# # otherwise it get rids of cluster 28..32
# Sycon <- Sycon %>%
#   MetacellsByGroups(group.by = c("orig.ident", "seurat_clusters"),
#                     min_cells = 10, k = 10,
#                     ident.group = 'seurat_clusters',
#                     slot = "scale.data",
#                     assay = "SCT",
#                     wgcna_name = "SCT_variableFeature")
# 
# # set up expr matrix
# Sycon <- Sycon %>%
#   SetDatExpr(assay = "SCT")
# 
# # test different soft power thresholds
# Sycon <- Sycon %>%
#   TestSoftPowers()
# plots <- Sycon %>%
#   PlotSoftPowers()
# print(patchwork::wrap_plots(plots, ncol=2))
# 
# # build the co-expr network
# # construct co-expression network
# Sycon <- Sycon %>%
#   ConstructNetwork(tom_outdir = "12_hdWGCNA/02_SCT_assay/",
#                    tom_name = "SCT_variableFeatures")
# 
# # compute module eigengenes and connectivity
# Sycon <- Sycon %>%
#   ModuleEigengenes()
# Sycon <- Sycon %>%
#   ModuleConnectivity()
# 
# # make a featureplot of hMEs for each module
# # here, plotted ranges are symmetrical (as hdWGCNA default settings)
# SCT_featurePlot_list <- Sycon %>%
#   ModuleFeaturePlot(features = "hMEs",
#                     order = TRUE)
# 
# print(patchwork::wrap_plots(SCT_featurePlot_list, ncol = 4))
# 
# MEs <- Sycon %>%
#   GetMEs(harmonized = TRUE)
# 
# modules <- GetModules(Sycon)
# mods <- levels(modules$module); mods <- mods[mods != 'grey']
# 
# modules <- GetModules(Sycon)
# mods <- levels(modules$module); mods <- mods[mods != 'grey']
# 
# Sycon@meta.data <- cbind(Sycon@meta.data, MEs)
# 
# 
# plot_module_UMAP("blue")
# 
# umap_list <- list()
# for (i in levels(Sycon@misc$SCT$wgcna_degrees$module)) {
#   if (i != "grey") {
#     umap_list[[i]] <- plot_module_UMAP(i) 
#   }
# }
# 
# umaps <- print(patchwork::wrap_plots(umap_list, ncol = 4))
# 
# ggsave("12_hdWGCNA/gene_module_umaps.pdf",
#        umaps,  device = cairo_pdf,
#        dpi = 300, height = 6, width = 12, units = ("in"), bg = 'white'
# )
# 
# ggsave("12_hdWGCNA/gene_module_umaps.png",
#        umaps,  device = "png",
#        dpi = 300, height = 6, width = 12, units = ("in"), bg = 'white'
# )
# 
# dotplot <- Sycon %>%
#   DotPlot(features = mods, group.by = 'seurat_clusters')
# # RotatedAxis()
# 
# dotplot_custom <- dotplot$data %>%
#   mutate(id = factor(id,
#                      levels = c("32", "27", "24", "19", "17",
#                                 "30", "25", "15",
#                                 "26", "21",
#                                 "23", "18", "16", "14", "11",
#                                 "31", "29", "28", "22", "20", "13", "12",
#                                 seq(10, 0, -1))),
#          features.plot = factor(features.plot,
#                                 levels = c("blue", "turquoise", "brown", "green", "yellow",
#                                            "pink", "red", "black"))) %>%
#   filter(pct.exp > 0) %>%
#   ggplot(aes(features.plot, id)) +
#   geom_point(aes(size = pct.exp, fill = avg.exp.scaled), shape = 21) +
#   
#   scale_y_discrete(limits = rev) +
#   scale_size_continuous(range = c(1,7)) +
#   
#   xlab("gene co-expression\nmodules") +
#   ylab("Cell clusters") +
#   labs(size = "Percent\nexpressed", fill = "Average\nexpression") +
#   
#   scale_fill_viridis_c(option = "plasma") +
#   
#   theme(strip.placement = "outside",
#         strip.background = element_blank(),
#         strip.text = element_text(face = "bold", size = 8),
#         strip.clip = "off",
#         axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
#         panel.border = element_rect(color = "#4f4f4f", fill = NA, linewidth = 0.8),
#         panel.background = element_blank(),
#         panel.grid = element_line(color = "#dbdbdb"))
# 
# dotplot_custom
# 
# ggsave("12_hdWGCNA/gene_module_dotplot.pdf",
#        dotplot_custom,  device = cairo_pdf,
#        dpi = 300, height = 6, width = 12, units = ("in"), bg = 'white'
# )
# 
# ggsave("12_hdWGCNA/gene_module_dotplot.png",
#        dotplot_custom,  device = "png",
#        dpi = 300, height = 6, width = 12, units = ("in"), bg = 'white'
# )