#!/usr/bin/env Rscript

setwd("/data/evassvis/fn76/sycon_clusterAnnotation")

suppressMessages({library(tidyverse)
  library(hdWGCNA)
  library(SeuratExtend)})


######################
#     LOAD INPUT     #
######################

load("00_input/Sycon_Seuratv4.Rdata")

Sycon <- UpdateSeuratObject(Sycon)


##########################################
#     hdWGCNA with SCTransformed data    #
##########################################

# THE DEVELOPERS DO NOT SUGGEST USING SCTransformed DATA

Sycon <- Sycon %>%
  SetupForWGCNA(features = VariableFeatures(Sycon),
                wgcna_name = "SCT")

Sycon <- Sycon %>%
  MetacellsByGroups(group.by = c("orig.ident", "seurat_clusters"),
                    reduction = 'harmony', # select the dimensionality reduction to perform KNN on
                    min_cells = 50,
                    ident.group = 'seurat_clusters', # set the Idents of the metacell seurat object
                    slot = "scale.data",
                    assay = "SCT"
  )

# set up expr matrix
Sycon <- Sycon %>%
  SetDatExpr(# group_name = "32", # the name of the group of interest in the group.by column
    # use_metacells = TRUE,
    # group.by = "seurat_clusters", # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
    assay = 'SCT' # using RNA assay
    # layer = 'counts' # using normalized data
  )

# test different soft power thresholds
Sycon <- Sycon %>%
  TestSoftPowers()
plots <- Sycon %>%
  PlotSoftPowers()
print(patchwork::wrap_plots(plots, ncol=2))

# build the co-expr network
Sycon <- Sycon %>%
  ConstructNetwork(tom_name = "SCT_cells",
                   overwrite_tom = TRUE
  )

# compute module eigengenes and connectivity
Sycon <- Sycon %>%
  ModuleEigengenes()
Sycon <- Sycon %>%
  ModuleConnectivity()

PlotDendrogram(Sycon, main = 'Sycon hdWGCNA Dendrogram')

PlotKMEs(Sycon, ncol=5)

module_features_plots <- Sycon %>%
  ModuleFeaturePlot(features = 'hMEs', # plot the hMEs
                    order = TRUE # order so the points with highest hMEs are on top
  )
print(patchwork::wrap_plots(module_features_plots, ncol = 4))

MEs <- Sycon %>%
  GetMEs(harmonized = TRUE)

modules <- GetModules(Sycon)
mods <- levels(modules$module); mods <- mods[mods != 'grey']

modules <- GetModules(Sycon)
mods <- levels(modules$module); mods <- mods[mods != 'grey']

Sycon@meta.data <- cbind(Sycon@meta.data, MEs)

plot_module_UMAP <- function(module_name) {
  
  umap <- Sycon@reductions$umap@cell.embeddings %>%
    as_tibble(rownames = NA) %>%
    rownames_to_column() %>%
    left_join(Sycon@meta.data %>% rownames_to_column()) %>%
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

DefaultAssay(Sycon) <- "RNA"

sycon_standard <- Sycon %>%
  DietSeurat(assays = "RNA", dimreducs = NULL)

colnames(sycon_standard@meta.data)[7] <- "seurat_clusters_og"

# standard seurat pipeline
sycon_standard <- NormalizeData(sycon_standard)
sycon_standard <- FindVariableFeatures(sycon_standard)
sycon_standard <- ScaleData(sycon_standard)
sycon_standard <- RunPCA(sycon_standard)

ElbowPlot(sycon_standard)

sycon_standard <- sycon_standard %>%
  FindNeighbors(dims = 1:15, reduction = "pca") %>%
  FindClusters()

sycon_standard[["RNA"]] <- split(sycon_standard[["RNA"]],
                                 f = sycon_standard$orig.ident)

sycon_standard <- IntegrateLayers(object = sycon_standard, method = HarmonyIntegration,
                                  orig.reduction = "pca", new.reduction = "harmony.int") #shit

sycon_standard <- IntegrateLayers(object = sycon_standard, method = CCAIntegration,
                                  orig.reduction = "pca", new.reduction = "cca.int")

sycon_standard <- IntegrateLayers(object = sycon_standard, method = RPCAIntegration,
                                  orig.reduction = "pca", new.reduction = "rpca.int")

sycon_standard[["RNA"]] <- JoinLayers(sycon_standard[["RNA"]])

sycon_standard <- sycon_standard %>%
  FindNeighbors(dims = 1:15, reduction = "rpca.int") %>%
  FindClusters()

sycon_standard <- sycon_standard %>%
  RunUMAP(dims = 1:15, reduction = "rpca.int")

DimPlot(sycon_standard,
        group.by = "orig.ident"
)

# standard hdWGCNA pipeline
sycon_standard <- SetupForWGCNA(
  sycon_standard,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "wgcna_standard" # the name of the hdWGCNA experiment
)

sycon_standard <- MetacellsByGroups(
  seurat_obj = sycon_standard,
  group.by = c("seurat_clusters_og", "orig.ident"), # specify the columns in seurat_obj@meta.data to group by
  reduction = "rpca.int", # select the dimensionality reduction to perform KNN on
  min_cells = 50,
  ident.group = "seurat_clusters_og" # set the Idents of the metacell seurat object
)

# sycon_metacell <- GetMetacellObject(sycon_standard)
# 
# sycon_metacell <- NormalizeData(sycon_metacell)
# sycon_metacell <- FindVariableFeatures(sycon_metacell)
# sycon_metacell <- ScaleData(sycon_metacell)
# sycon_metacell <- RunPCA(sycon_metacell)
# sycon_standard <- RunHarmonyMetacells(sycon_standard, group.by.vars = 'orig.ident')
# sycon_standard <- RunUMAPMetacells(sycon_standard, reduction = 'harmony', dims = 1:15)
# 
# DimPlotMetacells(sycon_standard, group.by = 'seurat_clusters_og') |
#   DimPlotMetacells(sycon_standard, group.by = 'orig.ident')

sycon_standard <- NormalizeMetacells(sycon_standard)
sycon_standard <- ScaleMetacells(sycon_standard,
                                 features = VariableFeatures(sycon_standard))
sycon_standard <- RunPCAMetacells(sycon_standard,
                                  features = VariableFeatures(sycon_standard))
sycon_standard <- RunHarmonyMetacells(sycon_standard,
                                      group.by.vars = "orig.ident")
sycon_standard <- RunUMAPMetacells(sycon_standard,
                                   reduction = 'harmony', dims = 1:15)

sycon_standard <- SetDatExpr(
  sycon_standard,
  # group_name = "INH", # the name of the group of interest in the group.by column
  # group.by = 'seurat_clusters_og', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  layer = 'data' # using normalized data
)
# Test different soft powers:
sycon_standard <- TestSoftPowers(
  sycon_standard,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plotSoftPower_list <- PlotSoftPowers(sycon_standard)

# assemble with patchwork
patchwork::wrap_plots(plotSoftPower_list, ncol = 2)

# construct co-expression network:
sycon_standard <- ConstructNetwork(
  sycon_standard
  # tom_name = 'INH' # name of the topoligical overlap matrix written to disk
)

# compute all MEs in the full single-cell dataset
sycon_standard <- ModuleEigengenes(
  sycon_standard,
  group.by.vars = "orig.ident"
)

# compute eigengene-based connectivity (kME):
sycon_standard <- ModuleConnectivity(
  sycon_standard,
  group.by = "seurat_clusters_og"
  # group_name = 'INH'
)

# make a featureplot of hMEs for each module
featurePlot_list <- ModuleFeaturePlot(
  sycon_standard,
  features = "hMEs", # plot the hMEs
  order = TRUE # order so the points with highest hMEs are on top
)

# stitch together with patchwork
patchwork::wrap_plots(featurePlot_list,
                      ncol = 3)

gene_list <- Sycon@misc$SCT$wgcna_modules %>%
  group_by(module) %>%
  summarise(genes = list(gene_name), .groups = "drop") %>%
  { setNames(.$genes, .$module)}

for (module in names(gene_list)) {
  file_path <- file.path("12_hdWGCNA", paste0(module, "_module_genes.ls"))
  writeLines(gene_list[[module]],
             file_path)
}

