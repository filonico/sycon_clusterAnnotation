#!/usr/bin/env Rscript

setwd("/data/evassvis/fn76/sycon_clusterAnnotation")

library(Seurat)
library(reticulate)
library(scCustomize)
library(tidyverse)
library(metacell)
# library(sceasy)

# initialize python from the conda env
reticulate::use_condaenv("RSeurat_env")

# import the anndata module to converto to/from .h5ad
ad <- reticulate::import("anndata", convert = FALSE)


#####################
#     I/O FILES     #
#####################

output_path <- "./04_preprocessed_scRNAseqs/"
input_path <- "./00_input/scRNAseqs_rawTables/"

# Sycon ciliatum files
scil_files <- list(Rdata = "00_input/Sycon_Seuratv4.Rdata",
                   out_h5ad = "Scil_cellFiltered.h5ad")

# Amphimedon queenslandica files
aque_files <- list(UMItable = paste0(input_path, "A_queenslandica/GSM3021561_Aque_Amphimedon_queenslandica_adult_UMI_table.txt.gz"),
                   metacell_assignments = paste0(input_path, "A_queenslandica/Amphimedon_adult_metacell_assignments.tsv"),
                   metacell_annotation = paste0(input_path, "A_queenslandica/Aque_adult_metacell_annotation.tsv"),
                   out_h5ad = "Aque_cellFiltered.h5ad")

# Spongilla lacustris files
slac_files <- list(UMItable = paste0(input_path, "S_lacustris/GSE134912_Slac_spongilla_10x_count_matrix_edited.txt.gz"),
                   cell_metadata = paste0(input_path, "S_lacustris/spongilla_cell_metadata_newNames.tsv"),
                   out_h5ad = "Slac_cellFiltered.h5ad")

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


#######################
#     S. CILIATUM     #
#######################

load(scil_files$Rdata)

scil_S.object <- Seurat::CreateSeuratObject(counts = Seurat::GetAssayData(object = Sycon,
                                                                          assay = "RNA",
                                                                          layer = 'counts'),
                                            project = "scil", min.cells = 0, min.features = 0) %>%
  Seurat::AddMetaData(metadata = Sycon@meta.data)

plot_featuresVScounts(scil_S.object)

# get batches
S231 <- scil_S.object %>%
  subset(orig.ident == "S231")
S232 <- scil_S.object %>%
  subset(orig.ident == "S232")
S251 <- scil_S.object %>%
  subset(orig.ident == "S251")
S252 <- scil_S.object %>%
  subset(orig.ident == "S252")

S231_S232 <- scil_S.object %>%
  subset(orig.ident %in% c("S231","S232"))
S251_S252 <- scil_S.object %>%
  subset(orig.ident %in% c("S251","S252"))

# save umi matrix
S231@assays$RNA$counts %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "GeneID") %>%
  write.table(file = "08_bonsai/00_umis_metadata/Scil_S231_raw_UMImatrix.tsv",
              col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
S232@assays$RNA$counts %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "GeneID") %>%
  write.table(file = "08_bonsai/00_umis_metadata/Scil_S232_raw_UMImatrix.tsv",
              col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
S251@assays$RNA$counts %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "GeneID") %>%
  write.table(file = "08_bonsai/00_umis_metadata/Scil_S251_raw_UMImatrix.tsv",
              col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
S252@assays$RNA$counts %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "GeneID") %>%
  write.table(file = "08_bonsai/00_umis_metadata/Scil_S252_raw_UMImatrix.tsv",
              col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

S231_S232@assays$RNA$counts %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "GeneID") %>%
  write.table(file = "08_bonsai/00_umis_metadata/Scil_S231_S232_raw_UMImatrix.tsv",
              col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

S251_S252@assays$RNA$counts %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "GeneID") %>%
  write.table(file = "08_bonsai/00_umis_metadata/Scil_S251_S252_raw_UMImatrix.tsv",
              col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

# save SCT umi matrix
Sycon@assays$SCT@counts %>%
  as_tibble(rownames = NA) %>%
  # rownames_to_column(var = "cellID") %>%
  write.table(file = "08_bonsai/00_umis_metadata/Scil_SCT_UMImatrix.tsv",
              col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

Sycon@assays$SCT@counts %>%
  as_tibble(rownames = NA) %>%
  rownames() %>%
  write.table(file = "08_bonsai/00_umis_metadata/Scil_SCT_geneID.tsv",
              col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

Sycon@assays$SCT@counts %>%
  as_tibble(rownames = NA) %>%
  colnames() %>%
  write.table(file = "08_bonsai/00_umis_metadata/Scil_SCT_cellID.tsv",
              col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")


# save metadata
Sycon@meta.data %>%
  rownames_to_column(var = "cellID") %>%
  write.table(file = "08_bonsai/00_umis_metadata/Scil_metadata.tsv",
              col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")



############################
#     A. QUEENSLANDICA     #
############################

# create the Seurat object from UMI table
aque_S.object <- read.table(file = gzfile(aque_files$UMItable),
                            header = TRUE, sep = "\t") %>%
  Seurat::CreateSeuratObject(project = "aque_adult", min.cells = 0, min.features = 0,
                             meta.data = read.table(aque_files$metacell_assignments, header = TRUE, sep = "\t") %>%
                               left_join(read.table(aque_files$metacell_annotation, header = TRUE, sep = "\t")) %>%
                               rename(cell_type = ID) %>%
                               column_to_rownames(var = "Cell"))

# plot features and counts before cell filtering
aque_beforeFiltering <- plot_featuresVScounts(aque_S.object)

# plot features and counts after filtering
aque_afterFiltering <- plot_featuresVScounts(aque_S.object,
                                             nCount_min = 200, nCount_max = 11000,
                                             nFeature_min = 200, nFeature_max = 4000)

# save umi matrix
aque_afterFiltering$S.object@assays$RNA$counts %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "GeneID") %>%
  write.table(file = "08_bonsai/00_umis_metadata/Aque_raw_UMImatrix.tsv",
              col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

# save metadata
aque_afterFiltering$S.object@meta.data %>%
  rownames_to_column(var = "cellID") %>%
  write.table(file = "08_bonsai/00_umis_metadata/Aque_metadata.tsv",
              col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")


########################
#     S. LACUSTRIS     #
########################

# create the Seurat object
slac_S.object <- read.table(file = gzfile(slac_files$UMItable),
                            header = TRUE, sep = "\t") %>%
  Seurat::CreateSeuratObject(project = "slac", min.cells = 0, min.features = 0,
                             
                             # add meta data with pre-computed cell clusters and annotations
                             meta.data = read.table(slac_files$cell_metadata, header = TRUE, sep = "\t") %>%
                               select(cell, clusterID, cell_type, cell_type_newName) %>%
                               mutate(cell = str_replace_all(cell, c('-' = '.'))) %>%
                               column_to_rownames(var = "cell"))

# plot features and counts before cell filtering
slac_beforeFiltering <- plot_featuresVScounts(slac_S.object)

# plot features and counts after cell filtering
slac_afterFiltering <- plot_featuresVScounts(slac_S.object,
                                             nCount_min = 200, nCount_max = 40000,
                                             nFeature_min = 200, nFeature_max = 6000)

# save UMI matrix
slac_afterFiltering$S.object@assays$RNA$counts %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "GeneID") %>%
  write.table(file = "08_bonsai/00_umis_metadata/Slac_raw_UMImatrix.tsv",
              col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

# save metadata
slac_afterFiltering$S.object@meta.data %>%
  rownames_to_column(var = "cellID") %>%
  write.table(file = "08_bonsai/00_umis_metadata/Slac_metadata.tsv",
              col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")