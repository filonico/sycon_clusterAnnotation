#!/usr/bin/env Rscript

setwd("/data/evassvis/fn76/sycon/sycon_clusterAnnotation")

suppressMessages({library(tidyverse)
  library(hdWGCNA)
  library(SeuratExtend)})


######################
#     LOAD INPUT     #
######################

Sycon_blobOnly <- readRDS("13_recluster_blob/Sycon_blobOnly.Rds") %>%
  UpdateSeuratObject()


##########################################
#     hdWGCNA with SCTransformed data    #
##########################################

# THE DEVELOPERS DO NOT SUGGEST USING SCTransformed DATA

Sycon_blobOnly <- Sycon_blobOnly %>%
  SetupForWGCNA(features = VariableFeatures(Sycon_blobOnly),
                wgcna_name = "SCT")

Sycon_blobOnly <- Sycon_blobOnly %>%
  MetacellsByGroups(group.by = "seurat_clusters_original",
                    #group.by = c("orig.ident", "seurat_clusters_original"),
                    # reduction = 'harmony', # select the dimensionality reduction to perform KNN on
                    min_cells = 50,
                    ident.group = 'seurat_clusters_original', # set the Idents of the metacell seurat object
                    slot = "scale.data",
                    assay = "SCT"
  )

# set up expr matrix
Sycon_blobOnly <- Sycon_blobOnly %>%
  SetDatExpr(# group_name = "32", # the name of the group of interest in the group.by column
    # use_metacells = TRUE,
    # group.by = "seurat_clusters", # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
    assay = 'SCT' # using RNA assay
    # layer = 'counts' # using normalized data
  )

# test different soft power thresholds
Sycon_blobOnly <- Sycon_blobOnly %>%
  TestSoftPowers()
plots <- Sycon_blobOnly %>%
  PlotSoftPowers()
print(patchwork::wrap_plots(plots, ncol=2))

# build the co-expr network
Sycon_blobOnly <- Sycon_blobOnly %>%
  ConstructNetwork(tom_name = "SCT_cells",
                   overwrite_tom = TRUE
  )

# compute module eigengenes and connectivity
Sycon_blobOnly <- Sycon_blobOnly %>%
  ModuleEigengenes()
Sycon_blobOnly <- Sycon_blobOnly %>%
  ModuleConnectivity()

PlotDendrogram(Sycon_blobOnly, main = 'Sycon hdWGCNA Dendrogram')

PlotKMEs(Sycon_blobOnly, ncol=5)

module_features_plots <- Sycon_blobOnly %>%
  ModuleFeaturePlot(features = 'hMEs', # plot the hMEs
                    order = TRUE # order so the points with highest hMEs are on top
  )
print(patchwork::wrap_plots(module_features_plots, ncol = 2))

MEs <- Sycon %>%
  GetMEs(harmonized = TRUE)

modules <- GetModules(Sycon)
mods <- levels(modules$module); mods <- mods[mods != 'grey']

modules <- GetModules(Sycon)
mods <- levels(modules$module); mods <- mods[mods != 'grey']

Sycon@meta.data <- cbind(Sycon@meta.data, MEs)

plot_module_UMAP <- function(module_name) {
  
  umap <- Sycon_blobOnly@reductions$umap@cell.embeddings %>%
    as_tibble(rownames = NA) %>%
    rownames_to_column() %>%
    left_join(Sycon_blobOnly@meta.data %>% rownames_to_column()) %>%
    select(rowname, UMAP_1, UMAP_2, all_of(module_name)) %>%
    pivot_longer(-c(rowname, UMAP_1, UMAP_2), names_to = "module", values_to = "expr") %>%
    arrange(expr) %>%
    # group_by(module) %>%
    
    ggplot(aes(UMAP_1, UMAP_2)) +
    # geom_point(col = "grey90", size = 0.7) +
    geom_point(aes(col = expr), #col = module_name,
               # alpha = 0.5, 
               size = 0.7
               ) +
    
    scale_color_gradient2(low="grey75", mid="grey95", high = module_name,
                          midpoint = 0) +
    # scale_alpha_continuous(range = c(0, 1)) +
    
    ggtitle(module_name) +
    
    # facet_wrap(~module) +
    # theme_umap_arrows() +
    coord_fixed() +
    
    theme(axis.text = element_text(size = 7),
          axis.title = element_text(size = 8),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          panel.border = element_rect(color = "#4f4f4f", fill = NA, linewidth = 0.8),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          legend.position = "none"
          )
  
  return(umap)
}

plot_module_UMAP("blue")

umap_list <- list()
for (i in levels(Sycon@misc$SCT$wgcna_degrees$module)) {
  if (i != "grey") {
    umap_list[[i]] <- plot_module_UMAP(i) 
  }
}

umaps <- print(patchwork::wrap_plots(umap_list, ncol = 4))

ggsave("12_hdWGCNA/gene_module_umaps.pdf",
       umaps,  device = cairo_pdf,
       dpi = 300, height = 6, width = 12, units = ("in"), bg = 'white'
)

ggsave("12_hdWGCNA/gene_module_umaps.png",
       umaps,  device = "png",
       dpi = 300, height = 6, width = 12, units = ("in"), bg = 'white'
)

dotplot <- Sycon %>%
  DotPlot(features = mods, group.by = 'seurat_clusters')
# RotatedAxis()

dotplot_custom <- dotplot$data %>%
  mutate(id = factor(id,
                     levels = c("32", "27", "24", "19", "17",
                                "30", "25", "15",
                                "26", "21",
                                "23", "18", "16", "14", "11",
                                "31", "29", "28", "22", "20", "13", "12",
                                seq(10, 0, -1))),
         features.plot = factor(features.plot,
                                levels = c("blue", "turquoise", "brown", "green", "yellow",
                                           "pink", "red", "black"))) %>%
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

dotplot_custom

ggsave("12_hdWGCNA/gene_module_dotplot.pdf",
       dotplot_custom,  device = cairo_pdf,
       dpi = 300, height = 6, width = 12, units = ("in"), bg = 'white'
)

ggsave("12_hdWGCNA/gene_module_dotplot.png",
       dotplot_custom,  device = "png",
       dpi = 300, height = 6, width = 12, units = ("in"), bg = 'white'
)

Sycon <- Sycon %>%
  RunModuleUMAP(n_hubs = 10, # number of hub genes to include for the UMAP embedding
                n_neighbors = 15, # neighbors parameter for UMAP
                min_dist = 0.1 # min distance between points in UMAP space
  )

# get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(Sycon)

# plot with ggplot
ggplot(umap_df, aes(UMAP1, UMAP2)) +
  geom_point(color=umap_df$color, size=umap_df$kME*2) +
  umap_theme()

options(future.globals.maxSize = 8000 * 1024^2)

Sycon %>%
  ModuleUMAPPlot(edge.alpha = 0.25,
                 sample_edges = TRUE,
                 edge_prop = 0.1, # proportion of edges to sample (20% here)
                 label_hubs = 2, # how many hub genes to plot per module?
                 keep_grey_edges = FALSE
  )


##########################################
#     hdWGCNA with standard pipeline     #
##########################################

DefaultAssay(Sycon_blobOnly) <- "RNA"

sycon_standardWGCNAworkflow <- Sycon_blobOnly %>%
  UpdateSeuratObject() %>%
  DietSeurat(assays = "RNA", dimreducs = NULL)

# standard seurat pipeline
sycon_standardWGCNAworkflow <- NormalizeData(sycon_standardWGCNAworkflow)
sycon_standardWGCNAworkflow <- FindVariableFeatures(sycon_standardWGCNAworkflow)
sycon_standardWGCNAworkflow <- ScaleData(sycon_standardWGCNAworkflow)
sycon_standardWGCNAworkflow <- RunPCA(sycon_standardWGCNAworkflow)

ElbowPlot(sycon_standardWGCNAworkflow)

DimPlot(sycon_standardWGCNAworkflow,
        group.by = "orig.ident")

sycon_standardWGCNAworkflow <- sycon_standardWGCNAworkflow %>%
  FindNeighbors(dims = 1:15, reduction = "pca") %>%
  FindClusters()

sycon_standardWGCNAworkflow <- sycon_standardWGCNAworkflow %>%
  RunUMAP(dims = 1:15)

DimPlot(sycon_standardWGCNAworkflow,
        group.by = "orig.ident")

sycon_standardWGCNAworkflow[["RNA"]] <- split(sycon_standardWGCNAworkflow[["RNA"]],
                                              f = sycon_standardWGCNAworkflow$orig.ident)

sycon_standardWGCNAworkflow <- IntegrateLayers(object = sycon_standardWGCNAworkflow, method = HarmonyIntegration,
                                  orig.reduction = "pca", new.reduction = "harmony.int") #shit

# this gives some error which I can't udnerstand
# sycon_standardWGCNAworkflow <- IntegrateLayers(object = sycon_standardWGCNAworkflow, method = CCAIntegration,
#                                   orig.reduction = "pca", new.reduction = "cca.int")

# rpca is way better
sycon_standardWGCNAworkflow <- IntegrateLayers(object = sycon_standardWGCNAworkflow, method = RPCAIntegration,
                                  orig.reduction = "pca", new.reduction = "rpca.int")

sycon_standardWGCNAworkflow[["RNA"]] <- JoinLayers(sycon_standardWGCNAworkflow[["RNA"]])

sycon_standardWGCNAworkflow <- sycon_standardWGCNAworkflow %>%
  FindNeighbors(dims = 1:15, reduction = "rpca.int") %>%
  FindClusters() %>%
  RunUMAP(dims = 1:15,
          reduction = "rpca.int")

DimPlot(sycon_standardWGCNAworkflow,
        group.by = "orig.ident",
        reduction = "umap") # rpca is way better

# standard hdWGCNA pipeline
sycon_standardWGCNAworkflow <- SetupForWGCNA(
  sycon_standardWGCNAworkflow,
  # gene_select = "fraction", # the gene selection approach
  # fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "wgcna_standard" # the name of the hdWGCNA experiment
)

sycon_standardWGCNAworkflow <- MetacellsByGroups(
  seurat_obj = sycon_standardWGCNAworkflow,
  # group.by = c("seurat_clusters", "orig.ident"), # specify the columns in seurat_obj@meta.data to group by
  reduction = "rpca.int", # select the dimensionality reduction to perform KNN on
  min_cells = 50,
  ident.group = "seurat_clusters" # set the Idents of the metacell seurat object
)

sycon_standardWGCNAworkflow <- NormalizeMetacells(sycon_standardWGCNAworkflow)
sycon_standardWGCNAworkflow <- ScaleMetacells(sycon_standardWGCNAworkflow,
                                              features = VariableFeatures(sycon_standardWGCNAworkflow))
sycon_standardWGCNAworkflow <- RunPCAMetacells(sycon_standardWGCNAworkflow,
                                               features = VariableFeatures(sycon_standardWGCNAworkflow))
# sycon_standardWGCNAworkflow <- RunHarmonyMetacells(sycon_standardWGCNAworkflow,
#                                                    group.by.vars = "orig.ident")
sycon_standardWGCNAworkflow <- RunUMAPMetacells(sycon_standardWGCNAworkflow,
                                                # reduction = 'harmony',
                                                dims = 1:15)

sycon_standardWGCNAworkflow <- SetDatExpr(
  sycon_standardWGCNAworkflow,
  # group_name = "INH", # the name of the group of interest in the group.by column
  # group.by = 'seurat_clusters_og', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  layer = 'data' # using normalized data
)

# Test different soft powers:
sycon_standardWGCNAworkflow <- TestSoftPowers(
  sycon_standardWGCNAworkflow
  # networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plotSoftPower_list <- PlotSoftPowers(sycon_standardWGCNAworkflow)

# assemble with patchwork
patchwork::wrap_plots(plotSoftPower_list, ncol = 2)

# construct co-expression network:
sycon_standardWGCNAworkflow <- ConstructNetwork(
  sycon_standardWGCNAworkflow,
  tom_name = 'TOM_blobOnly', # name of the topoligical overlap matrix written to disk
  overwrite_tom = TRUE
)

# compute all MEs in the full single-cell dataset
sycon_standardWGCNAworkflow <- ModuleEigengenes(
  sycon_standardWGCNAworkflow
  # group.by.vars = "orig.ident"
)

# compute eigengene-based connectivity (kME):
sycon_standardWGCNAworkflow <- ModuleConnectivity(
  sycon_standardWGCNAworkflow,
  # group.by = "seurat_clusters",
  group_name = 'TOM_blobOnly'
)

# make a featureplot of hMEs for each module
featurePlot_list <- ModuleFeaturePlot(
  sycon_standardWGCNAworkflow,
  features = "hMEs", # plot the hMEs
  order = TRUE # order so the points with highest hMEs are on top
)

# stitch together with patchwork
umaps <- patchwork::wrap_plots(featurePlot_list,
                      ncol = 2)

ggsave("13_recluster_blob/03_hdWGCNA/gene_module_umaps.pdf",
       umaps,  device = cairo_pdf,
       dpi = 300, height = 6, width = 12, units = ("in"), bg = 'white'
)

ggsave("13_recluster_blob/03_hdWGCNA/gene_module_umaps.png",
       umaps,  device = "png",
       dpi = 300, height = 6, width = 12, units = ("in"), bg = 'white'
)

gene_list <- Sycon@misc$SCT$wgcna_modules %>%
  group_by(module) %>%
  summarise(genes = list(gene_name), .groups = "drop") %>%
  { setNames(.$genes, .$module)}

for (module in names(gene_list)) {
  file_path <- file.path("12_hdWGCNA", paste0(module, "_module_genes.ls"))
  writeLines(gene_list[[module]],
             file_path)
}

