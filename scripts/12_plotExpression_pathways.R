setwd("/scratch/evassvis/fn76/ANALYSIS/origin_of_vision/")

library(Seurat)
library(SeuratExtend)
library(reticulate)
library(scCustomize)
library(tidyverse)
library(metacell)


#####################
#     I/O FILES     #
#####################

output_path <- "./04_preprocessed_scRNAseqs/"
input_path <- "./00_input/scRNAseqs_rawTables/"

scil_files <- list(Rdata = "00_input/Sycon_Seuratv4.Rdata",
                   out_h5ad = "Scil_cellFiltered.h5ad")

# Sycon ciliatum files
scil_files_reAnno <- list(in_h5ad = "00_input/Sycon.ReAnno.h5ad",
                          out_h5ad = "Scil_cellFiltered.h5ad")


#####################
#     FUNCTIONS     #
#####################

plot_featuresVScounts <- function(seuratObject,
                                  nCount_min, nCount_max,
                                  nFeature_min, nFeature_max) {
  
  # if none of the filtering parameter is set, then plot non-filtered data
  if (missing(nCount_min) && missing(nCount_max) && missing(nFeature_min) && missing(nFeature_max)) {
    seuratObject_toPlot <- seuratObject
  }
  
  # else, filter low-quality cell and then plot to data
  else {
    seuratObject_toPlot <- seuratObject %>%
      subset(subset = nFeature_RNA > nFeature_min & nFeature_RNA < nFeature_max &
               nCount_RNA > nCount_min & nCount_RNA < nCount_max)
  }
  
  # produce violin plot
  vlnplot <- VlnPlot(seuratObject_toPlot,
                     features = c("nFeature_RNA", "nCount_RNA"), group.by = "orig.ident")
  
  # produce scatter plot
  scatterplot <- FeatureScatter(seuratObject_toPlot,
                                feature1 = "nCount_RNA", feature2 = "nFeature_RNA") #+
  # scale_x_log10() +
  # scale_y_log10()
  
  # put the plots in a panel
  panel <- ggpubr::ggarrange(vlnplot, scatterplot)
  
  if (missing(nCount_min) && missing(nCount_max) && missing(nFeature_min) && missing(nFeature_max)) {
    return(panel)
  }
  else {
    return(list(S.object = seuratObject_toPlot, panel = panel))
  }
}


rename_cell_types <- function(seurat_metadata) {
  seurat_metadata <- seurat_metadata %>%
    mutate(cell_type_reduced = str_replace(cell_type, "_\\d+(\\.\\d+)?$", ""),
           cell_type_reduced = if_else(
             str_starts(cell_type_reduced, "peptidergic"),
             str_replace(cell_type_reduced, "_.*", ""),
             cell_type_reduced),
           cell_type_reduced = str_replace_all(cell_type_reduced, "_", " "))
  
  return(seurat_metadata)
}

#######################
#     S. CILIATUM     #
#######################

load(scil_files$Rdata)
scil_S.object <- Seurat::CreateSeuratObject(counts = Seurat::GetAssayData(object = Sycon,
                                                                          assay = "RNA",
                                                                          layer = 'counts'),
                                            project = "scil", min.cells = 0, min.features = 0) %>%
  Seurat::AddMetaData(metadata = Sycon@meta.data)


sceasy::convertFormat(scil_files_reAnno$in_h5ad, from = "anndata", to = "seurat", outFile = "00_input/Sycon.ReAnno.rds")
scil_S.object_reAnno <- readRDS("00_input/Sycon.ReAnno.rds")

# plot features and counts before cell filtering
scil_afterFiltering <- plot_featuresVScounts(scil_S.object_reAnno,
                                             0, 20000, 0, 3000)


####################################
#     DIMENSIONALITY REDUCTION     #
####################################

# perform data transformation and run PCA
scil_afterFiltering$S.object <- SCTransform(scil_afterFiltering$S.object) %>%
  SetIdent(value = 'N1.Int.Clusters') %>%
  RunPCA() #%>%

# run UMAP reduction
ElbowPlot(scil_afterFiltering$S.object, ndims = 50) # 20
scil_afterFiltering$S.object <- scil_afterFiltering$S.object %>%
  RunUMAP(dims = 1:20)

scil_umap_recoded <- scil_afterFiltering$S.object %>%
  SetIdent(value = 'N1.Int.Clusters') %>%
  DimPlot2(theme = theme_umap_arrows(), pt.size = 2, label = TRUE, repel = TRUE)

scil_umap_recoded

ggsave("06_variousPathway_genes/scil_umap_recodedCellTypes.pdf",
       scil_umap_recoded, device = cairo_pdf,
       dpi = 300, height = 7, width = 7, units = ("in"), bg = 'white')

scil_umap_og <- scil_afterFiltering$S.object %>%
  SetIdent(value = 'seurat_clusters') %>%
  DimPlot2(theme = theme_umap_arrows(), pt.size = 2, label = TRUE, repel = TRUE)

ggsave("06_variousPathway_genes/scil_umap_ogCellTypes.pdf",
       scil_umap_og, device = cairo_pdf,
       dpi = 300, height = 7, width = 7, units = ("in"), bg = 'white')

############################
#     FIND ALL MARKERS     #
############################

# find all markers
scil_allMarkers <- FindAllMarkers(scil_afterFiltering$S.object)

# write all markers to tables
scil_allMarkers %>%
  write.table(file = "07_allMarkers/scil_allMarkers.tsv", quote = FALSE, sep = '\t', row.names = FALSE)


#############################################
#     ANALYSE PHOTOTRANSDUCTION MARKERS     #
#############################################

scil_photoGenes <- read.table(file = "00_input/Scil_variousPathway_genes.tsv", header = FALSE, sep = "\t") %>%
  setNames(c("gene_description", "gene_ID", "pathway")) %>%
  filter(pathway == "visual transduction") %>%
  filter(!gene_ID %in% c("g2383", "g2865", "g7884", "g7972")) %>%
  group_by(gene_description) %>%
  summarize(gene_ID = list(gene_ID)) %>%
  with(setNames(gene_ID, gene_description))

# check expression per cell type
scil_dotplot_photoGenes_cellType <- scil_afterFiltering$S.object %>%
  SetIdent(value = 'seurat_clusters') %>%
  # SetIdent(value = 'N1.Int.Clusters') %>%
  DotPlot2(features = scil_photoGenes, flip = TRUE, )

scil_dotplot_photoGenes_cellType

# save dotplots
ggsave("06_variousPathway_genes/scil_photoGenes_dotplot_expression.pdf",
       plot = scil_dotplot_photoGenes_cellType, device = cairo_pdf,
       dpi = 300, height = 7, width = 12, units = ("in"), bg = 'white')

# set colors for expression
scil_list_photoGenes <- unlist(scil_photoGenes, use.names = FALSE)
scil_colorScale_photoGenes <- setNames(
  lapply(scil_list_photoGenes, function(x) c(low = "lightblue", high = "darkblue")),
  scil_list_photoGenes
)

# plot expression of all the monoaminergic markers
scil_photoGenes_expression <- scil_afterFiltering$S.object %>%
  DimPlot2(theme = NoAxes(),
           features = scil_list_photoGenes,
           cols = scil_colorScale_photoGenes,
           pt.size = 0.7, label = TRUE, repel = TRUE) +
  theme_umap_arrows() +
  patchwork::plot_annotation("Phototransduction-related genes")
scil_photoGenes_expression
  
# save panel
ggsave("06_variousPathway_genes/scil_umap_photoGenes.pdf",
       plot = scil_photoGenes_expression, device = cairo_pdf,
       dpi = 300, height = 12, width = 18, units = ("in"), bg = 'white')


#########################################
#     ANALYSE NOTCH PATHWAY MARKERS     #
#########################################

scil_notchGenes <- read.table(file = "00_input/Scil_variousPathway_genes.tsv", header = FALSE, sep = "\t") %>%
  setNames(c("gene_description", "gene_ID", "pathway")) %>%
  filter(pathway == "notch pathway") %>%
  filter(!gene_ID %in% c("g5692", "g2517")) %>%
  group_by(gene_description) %>%
  summarize(gene_ID = list(gene_ID)) %>%
  with(setNames(gene_ID, gene_description))

# check expression per cell type
scil_dotplot_notchGenes_cellType <- scil_afterFiltering$S.object %>%
  SetIdent(value = 'seurat_clusters') %>%
  # SetIdent(value = 'N1.Int.Clusters') %>%
  DotPlot2(features = scil_notchGenes, flip = TRUE, )

scil_dotplot_notchGenes_cellType

# save dotplots
ggsave("06_variousPathway_genes/scil_notchGenes_dotplot_expression.pdf",
       plot = scil_dotplot_notchGenes_cellType, device = cairo_pdf,
       dpi = 300, height = 7, width = 12, units = ("in"), bg = 'white')

# set colors for expression
scil_list_notchGenes <- unlist(scil_notchGenes, use.names = FALSE)
scil_colorScale_notchGenes <- setNames(
  lapply(scil_list_notchGenes, function(x) c(low = "lightblue", high = "darkblue")),
  scil_list_notchGenes
)

# plot expression of all the monoaminergic markers
scil_notchGenes_expression <- scil_afterFiltering$S.object %>%
  DimPlot2(theme = NoAxes(),
           features = scil_list_notchGenes,
           cols = scil_colorScale_notchGenes,
           pt.size = 0.7, label = TRUE, repel = TRUE) +
  theme_umap_arrows() +
  patchwork::plot_annotation("Notch-pathway genes")
scil_notchGenes_expression

# save panel
ggsave("06_variousPathway_genes/scil_umap_notchGenes.pdf",
       plot = scil_notchGenes_expression, device = cairo_pdf,
       dpi = 300, height = 12, width = 18, units = ("in"), bg = 'white')


#######################################
#     ANALYSE TGF PATHWAY MARKERS     #
#######################################

scil_tgfGenes <- read.table(file = "00_input/Scil_variousPathway_genes.tsv", header = FALSE, sep = "\t") %>%
  setNames(c("gene_description", "gene_ID", "pathway")) %>%
  filter(pathway == "TGF pathway") %>%
  # filter(!gene_ID %in% c("g5692", "g2517")) %>%
  group_by(gene_description) %>%
  summarize(gene_ID = list(gene_ID)) %>%
  with(setNames(gene_ID, gene_description))

# check expression per cell type
scil_dotplot_tgfGenes_cellType <- scil_afterFiltering$S.object %>%
  SetIdent(value = 'seurat_clusters') %>%
  # SetIdent(value = 'N1.Int.Clusters') %>%
  DotPlot2(features = scil_tgfGenes, flip = TRUE, )

scil_dotplot_tgfGenes_cellType

# save dotplots
ggsave("06_variousPathway_genes/scil_tgfGenes_dotplot_expression.pdf",
       plot = scil_dotplot_tgfGenes_cellType, device = cairo_pdf,
       dpi = 300, height = 7, width = 12, units = ("in"), bg = 'white')

# set colors for expression
scil_list_tgfGenes <- unlist(scil_tgfGenes, use.names = FALSE)
scil_colorScale_tgfGenes <- setNames(
  lapply(scil_list_tgfGenes, function(x) c(low = "lightblue", high = "darkblue")),
  scil_list_tgfGenes
)

# plot expression of all the monoaminergic markers
scil_tgfGenes_expression <- scil_afterFiltering$S.object %>%
  DimPlot2(theme = NoAxes(),
           features = scil_list_tgfGenes,
           cols = scil_colorScale_tgfGenes,
           pt.size = 0.7, label = TRUE, repel = TRUE) +
  theme_umap_arrows() +
  patchwork::plot_annotation("tgf-pathway genes")
scil_tgfGenes_expression

# save panel
ggsave("06_variousPathway_genes/scil_umap_tgfGenes.pdf",
       plot = scil_tgfGenes_expression, device = cairo_pdf,
       dpi = 300, height = 12, width = 18, units = ("in"), bg = 'white')


#######################################
#     ANALYSE WNT PATHWAY MARKERS     #
#######################################

scil_wntGenes <- read.table(file = "00_input/Scil_variousPathway_genes.tsv", header = FALSE, sep = "\t") %>%
  setNames(c("gene_description", "gene_ID", "pathway")) %>%
  filter(pathway == "Wnt pathway") %>%
  filter(!gene_ID %in% c("g13396", "g6493", "g13684", "g4399", "g778")) %>%
  group_by(gene_description) %>%
  summarize(gene_ID = list(gene_ID)) %>%
  with(setNames(gene_ID, gene_description))

# check expression per cell type
scil_dotplot_wntGenes_cellType <- scil_afterFiltering$S.object %>%
  SetIdent(value = 'seurat_clusters') %>%
  # SetIdent(value = 'N1.Int.Clusters') %>%
  DotPlot2(features = scil_wntGenes, flip = TRUE, )

scil_dotplot_wntGenes_cellType

# save dotplots
ggsave("06_variousPathway_genes/scil_wntGenes_dotplot_expression.pdf",
       plot = scil_dotplot_wntGenes_cellType, device = cairo_pdf,
       dpi = 300, height = 7, width = 12, units = ("in"), bg = 'white')

# set colors for expression
scil_list_wntGenes <- unlist(scil_wntGenes, use.names = FALSE)
scil_colorScale_wntGenes <- setNames(
  lapply(scil_list_wntGenes, function(x) c(low = "lightblue", high = "darkblue")),
  scil_list_wntGenes
)

# plot expression of all the monoaminergic markers
scil_wntGenes_expression <- scil_afterFiltering$S.object %>%
  DimPlot2(theme = NoAxes(),
           features = scil_list_wntGenes,
           cols = scil_colorScale_wntGenes,
           pt.size = 0.7, label = TRUE, repel = TRUE) +
  theme_umap_arrows() +
  patchwork::plot_annotation("Wnt-pathway genes")
scil_wntGenes_expression

# save panel
ggsave("06_variousPathway_genes/scil_umap_wntGenes.pdf",
       plot = scil_wntGenes_expression, device = cairo_pdf,
       dpi = 300, height = 12, width = 18, units = ("in"), bg = 'white')


###################################################
#     ANALYSE BASAL EPITHELIA PATHWAY MARKERS     #
###################################################

scil_basalEpitheliaGenes <- read.table(file = "00_input/Scil_variousPathway_genes.tsv", header = FALSE, sep = "\t") %>%
  setNames(c("gene_description", "gene_ID", "pathway")) %>%
  filter(pathway == "Basal epithelia pathway") %>%
  # filter(!gene_ID %in% c("g13396", "g6493", "g13684", "g4399", "g778")) %>%
  group_by(gene_description) %>%
  summarize(gene_ID = list(gene_ID)) %>%
  with(setNames(gene_ID, gene_description))

# check expression per cell type
scil_dotplot_basalEpitheliaGenes_cellType <- scil_afterFiltering$S.object %>%
  SetIdent(value = 'seurat_clusters') %>%
  # SetIdent(value = 'N1.Int.Clusters') %>%
  DotPlot2(features = scil_basalEpitheliaGenes, flip = TRUE, )

scil_dotplot_basalEpitheliaGenes_cellType

# save dotplots
ggsave("06_variousPathway_genes/scil_basalEpitheliaGenes_dotplot_expression.pdf",
       plot = scil_dotplot_basalEpitheliaGenes_cellType, device = cairo_pdf,
       dpi = 300, height = 7, width = 12, units = ("in"), bg = 'white')

# set colors for expression
scil_list_basalEpitheliaGenes <- unlist(scil_basalEpitheliaGenes, use.names = FALSE)
scil_colorScale_basalEpitheliaGenes <- setNames(
  lapply(scil_list_basalEpitheliaGenes, function(x) c(low = "lightblue", high = "darkblue")),
  scil_list_basalEpitheliaGenes
)

# plot expression of all the monoaminergic markers
scil_basalEpitheliaGenes_expression <- scil_afterFiltering$S.object %>%
  DimPlot2(theme = NoAxes(),
           features = scil_list_basalEpitheliaGenes,
           cols = scil_colorScale_basalEpitheliaGenes,
           pt.size = 0.7, label = TRUE, repel = TRUE) +
  theme_umap_arrows() +
  patchwork::plot_annotation("basalEpithelia-pathway genes")
scil_basalEpitheliaGenes_expression

# save panel
ggsave("06_variousPathway_genes/scil_umap_basalEpitheliaGenes.pdf",
       plot = scil_basalEpitheliaGenes_expression, device = cairo_pdf,
       dpi = 300, height = 12, width = 18, units = ("in"), bg = 'white')


#############################################
#     ANALYSE BIOMINERALIZATION MARKERS     #
#############################################

scil_biomineralGenes <- read.table(file = "00_input/Scil_variousPathway_genes.tsv", header = FALSE, sep = "\t") %>%
  setNames(c("gene_description", "gene_ID", "pathway")) %>%
  filter(pathway == "Biomineralization") %>%
  filter(!gene_ID %in% c("g782", "g11205")) %>%
  group_by(gene_description) %>%
  summarize(gene_ID = list(gene_ID)) %>%
  with(setNames(gene_ID, gene_description))

# check expression per cell type
scil_dotplot_biomineralGenes_cellType <- scil_afterFiltering$S.object %>%
  SetIdent(value = 'seurat_clusters') %>%
  # SetIdent(value = 'N1.Int.Clusters') %>%
  DotPlot2(features = scil_biomineralGenes, flip = TRUE, )

scil_dotplot_biomineralGenes_cellType

# save dotplots
ggsave("06_variousPathway_genes/scil_biomineralGenes_dotplot_expression.pdf",
       plot = scil_dotplot_biomineralGenes_cellType, device = cairo_pdf,
       dpi = 300, height = 7, width = 12, units = ("in"), bg = 'white')

# set colors for expression
scil_list_biomineralGenes <- unlist(scil_biomineralGenes, use.names = FALSE)
scil_colorScale_biomineralGenes <- setNames(
  lapply(scil_list_biomineralGenes, function(x) c(low = "lightblue", high = "darkblue")),
  scil_list_biomineralGenes
)

# plot expression of all the monoaminergic markers
scil_biomineralGenes_expression <- scil_afterFiltering$S.object %>%
  DimPlot2(theme = NoAxes(),
           features = scil_list_biomineralGenes,
           cols = scil_colorScale_biomineralGenes,
           pt.size = 0.7, label = TRUE, repel = TRUE) +
  theme_umap_arrows() +
  patchwork::plot_annotation("biomineral-pathway genes")
scil_biomineralGenes_expression

# save panel
ggsave("06_variousPathway_genes/scil_umap_biomineralGenes.pdf",
       plot = scil_biomineralGenes_expression, device = cairo_pdf,
       dpi = 300, height = 12, width = 18, units = ("in"), bg = 'white')
