#!/usr/bin/env Rscript

setwd("/data/evassvis/fn76/ANALYSIS/sycon_clusterAnnotation")

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

# scil_files_reAnno <- list(in_h5ad = "00_input/Sycon.ReAnno.h5ad",
#                           out_h5ad = "Scil_cellFiltered.h5ad")

# Amphimedon queenslandica files
aque_files <- list(UMItable = paste0(input_path, "A_queenslandica/GSM3021561_Aque_Amphimedon_queenslandica_adult_UMI_table.txt.gz"),
                   metacell_assignments = paste0(input_path, "A_queenslandica/Amphimedon_adult_metacell_assignments.tsv"),
                   metacell_annotation = paste0(input_path, "A_queenslandica/Aque_adult_metacell_annotation.tsv"),
                   out_h5ad = "Aque_cellFiltered.h5ad")

# # Clytia hemisphaerica files
# chem_CellRanger_dir <- "00_input/scRNAseqs_rawTables/C_hemisphaerica/"
# chem_h5ad_file <- "Chem_cellFiltered.h5ad"

# Hoilungia hongkongensis H13 files
hhon_files <- list(dir = paste0(input_path, "H_hongkongensis_H13/"),
                   Rda = "scdr_Hhon_it2",
                   mc = "scdr_Hhon_it4",
                   metacell_annotation = paste0(input_path, "H_hongkongensis_H13/annotation_mc.Hhon.it4.tsv"),
                   out_h5ad = "Hhon_cellFiltered.h5ad")

# Cladtertia collaboinventa (Hoilungia sp. H23) files
hh23_files <- list(dir = paste0(input_path, "H_sp_H23/"),
                   Rda = "scdr_HoiH23_it2",
                   mc = "scdr_HoiH23_it4",
                   metacell_annotation = paste0(input_path, "H_sp_H23/annotation_mc.HoiH23.it4.tsv"),
                   out_h5ad = "HH23_cellFiltered.h5ad")

# Hydra vulgaris files
hvul_files <- list(RDS = paste0(input_path, "H_vulgaris/Hvul_sc_UMI_counts.RDS"),
                   metacell_assignments = paste0(input_path, "H_vulgaris/Hvul_metacell_assignments.tsv"),
                   metacell_annotation = paste0(input_path, "H_vulgaris/Hvul_metacell_annotation.tsv"),
                   out_h5ad = "Hvul_cellFiltered.h5ad")

# Mnemiopsis leidyi files
mlei_files <- list(UMItable = paste0(input_path, "M_leidyi/GSM3021563_Mlei_Mnemiopsis_leidyi_UMI_table.txt.gz"),
                   metacell_assignments = paste0(input_path, "M_leidyi/Mnemiopsis_metacell_assignments.tsv"),
                   metacell_annotation = paste0(input_path, "M_leidyi/Mnemiopsis_metacell_annotation.tsv"),
                   out_h5ad = "Mlei_cellFiltered.h5ad")
                   
# Nematostella vectensis files
nvec_files <- list(RDS = paste0(input_path, "N_vectensis/Nvec_adult_sc_UMI_counts.RDS"),
                   metacell_assignments = paste0(input_path, "N_vectensis/Nematostella_adult_metacell_assignments.tsv"),
                   metacell_annotation = paste0(input_path, "N_vectensis/Nematostella_adult_metacell_annotation.tsv"),
                   out_h5ad = "Nvec_cellFiltered.h5ad")

# Spongilla lacustris files
slac_files <- list(UMItable = paste0(input_path, "S_lacustris/GSE134912_Slac_spongilla_10x_count_matrix_edited.txt.gz"),
                   cell_metadata = paste0(input_path, "S_lacustris/spongilla_cell_metadata_newNames.tsv"),
                   out_h5ad = "Slac_cellFiltered.h5ad")

# Stilophora pistillata files
spis_files <- list(RDS = paste0(input_path, "S_pistillata/Spis_adult_sc_UMI_counts.RDS"),
                   metacell_assignments = paste0(input_path, "S_pistillata/Spis_coral_metacell_assignments.tsv"),
                   metacell_annotation = paste0(input_path, "S_pistillata/Spis_coral_metacell_annotation.tsv"),
                   out_h5ad = "Spis_cellFiltered.h5ad")

# Trichoplax adhaerens files
tadh_files <- list(dir = paste0(input_path, "T_adhaerens/"),
                   Rda = "scdr_Tadh_it2",
                   mc = "scdr_Tadh_it4",
                   metacell_annotation = paste0(input_path, "T_adhaerens/annotation_mc.Tadh.it4.tsv"),
                   out_h5ad = "Tadh_cellFiltered.h5ad")

# Trichoplax sp. H2 files
trh2_files <- list(dir = paste0(input_path, "T_adhaerens_H2/"),
                   Rda = "scdr_TrH2_it2",
                   mc = "scdr_TrH2_it4",
                   metacell_annotation = paste0(input_path, "T_adhaerens_H2/annotation_mc.TrH2.it4.tsv"),
                   out_h5ad = "TrH2_cellFiltered.h5ad")

# Xenia sp. files
xesp_files <- list(RDS = paste0(input_path, "X_sp/Xesp_sc_UMI_counts.RDS"),
                   metacell_assignments = paste0(input_path, "X_sp/Xesp_metacell_assignments.tsv"),
                   metacell_annotation = paste0(input_path, "X_sp/Xesp_metacell_annotation.tsv"),
                   out_h5ad = "Xesp_cellFiltered.h5ad")


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

# function to convert from Seurat to AnnData object and save to file
fromSeurat_toAnndata <- function(s.object, outdir, outfile) {
  
  scCustomize::as.anndata(x = s.object, main_layer = "counts",
                          other_layers = NULL, file_path = outdir,
                          file_name = outfile)
}

# function to load placozoan scRNA-seq, all of which share the same file data structures
load_placozoa_data <- function(directory,
                               metacell_file,
                               Rda_file,
                               spID,
                               metacell_annotation) {
  
  # initialise metacell database
  metacell::scdb_init(directory, force_reinit = TRUE)
  
  # read in the metacell file
  mc <- metacell::scdb_mc(metacell_file)
  
  # read the .Rda data object and create a Seurat Object
  S.object <- metacell::scdb_mat(Rda_file)@mat %>%
    Seurat::CreateSeuratObject(project = spID, min.cells = 0, min.features = 0,
                               meta.data = enframe(mc@mc, value = "Cluster") %>%
                                 column_to_rownames(var = "name"))
  
  # add metacell annotation
  S.object@meta.data <- S.object@meta.data %>%
    rownames_to_column(var = "Cell") %>%
    left_join(read.table(metacell_annotation, header = TRUE, sep = "\t") %>%
                select(metacell, cell_type),
              by = join_by(Cluster == metacell)) %>%
    column_to_rownames(var = "Cell")
  
  return(S.object)
  
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

# sceasy::convertFormat(scil_files_reAnno$in_h5ad, from = "anndata", to = "seurat", outFile = "Sycon.ReAnno.rds")
# scil_S.object_reAnno <- readRDS("Sycon.ReAnno.rds")

# convert SeuratObject to anndata file
fromSeurat_toAnndata(scil_S.object,
                       output_path, scil_files$out_h5ad)

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

# convert SeuratObject to anndata file
fromSeurat_toAnndata(aque_afterFiltering$S.object,
                     output_path, aque_files$out_h5ad)


# ############################
# #     C. HEMISPHAERICA     #
# ############################
# 
# # NOT USING #
# 
# # create the Seurat object from CellRanger output
# chem_S.object <- Seurat::Read10X(chem_CellRanger_dir) %>%
#   Seurat::CreateSeuratObject(project = "chem", min.features = 0, min.cells = 0)
# 
# # plot features and counts before cell filtering
# chem_beforeFiltering <- plot_featuresVScounts(chem_S.object)
# 
# # plot features and counts before cell filtering
# chem_afterFiltering <- plot_featuresVScounts(chem_S.object,
#                                              nCount_min = 200, nCount_max = 2000,
#                                              nFeature_min = 200, nFeature_max = 750)
# 
# # convert SeuratObject to anndata file
# fromSeurat_toAnndata(chem_afterFiltering$S.object, output_path, chem_h5ad_file)


################################
#     H. HONGKONGENSIS H13     #
################################

# create the Seurat object fromt the UMI table
hhon_S.object <- load_placozoa_data(hhon_files$dir,
                                    hhon_files$mc,
                                    hhon_files$Rda,
                                    "hhon",
                                    hhon_files$metacell_annotation)

# plot features and counts before cell filtering
hhon_beforeFiltering <- plot_featuresVScounts(hhon_S.object)

# plot features and counts before cell filtering
hhon_afterFiltering <- plot_featuresVScounts(hhon_S.object,
                                             nCount_min = 200, nCount_max = 6500,
                                             nFeature_min = 200, nFeature_max = 1800)

# convert SeuratObject to anndata file
fromSeurat_toAnndata(hhon_afterFiltering$S.object,
                     output_path, hhon_files$out_h5ad)


######################
#     H. SP. H23     #
######################

# create the Seurat object fromt the UMI table
hh23_S.object <- load_placozoa_data(hh23_files$dir,
                                    hh23_files$mc,
                                    hh23_files$Rda,
                                    "hh23",
                                    hh23_files$metacell_annotation)

# plot features and counts before cell filtering
hh23_beforeFiltering <- plot_featuresVScounts(hh23_S.object)

# plot features and counts before cell filtering
hh23_afterFiltering <- plot_featuresVScounts(hh23_S.object,
                                             nCount_min = 200, nCount_max = 7500,
                                             nFeature_min = 200, nFeature_max = 2100)
# convert SeuratObject to anndata file
fromSeurat_toAnndata(hh23_afterFiltering$S.object,
                     output_path, hh23_files$out_h5ad)


########################
#      H. VULGARIS     #
########################

# create the Seurat object
hvul_S.object <- readRDS(hvul_files$RDS) %>%
  Seurat::CreateSeuratObject(project = "hvul", min.cells = 0, min.features = 0,
                             
                             # add meta data with pre-computed cell clusters and annotations
                             meta.data = read.table(hvul_files$metacell_assignments, header = TRUE, sep = "\t") %>%
                               left_join(read.table(hvul_files$metacell_annotation, header = TRUE, sep = "\t") %>%
                                           select(metacell, cell_type)) %>%
                               column_to_rownames(var = "cell"))

# plot features and counts before cell filtering
hvul_beforeFiltering <- plot_featuresVScounts(hvul_S.object)

# plot features and counts before cell filtering
hvul_afterFiltering <- plot_featuresVScounts(hvul_S.object,
                                             nCount_min = 200, nCount_max = 25000,
                                             nFeature_min = 200, nFeature_max = 4500)

# convert SeuratObject to anndata file
fromSeurat_toAnndata(hvul_S.object,
                     output_path, hvul_files$out_h5ad)


#####################
#     M. LEIDYI     #
#####################

# create the Seurat object
mlei_S.object <- read.table(file = gzfile(mlei_files$UMItable),
                            header = TRUE, sep = "\t") %>%
  Seurat::CreateSeuratObject(project = "mlei", min.cells = 0, min.features = 0,
                             
                             # add meta data with pre-computed cell clusters and annotations
                             meta.data = read.table(mlei_files$metacell_assignments, header = TRUE, sep = "\t") %>%
                               left_join(read.table(mlei_files$metacell_annotation, header = TRUE, sep = "\t"),
                                         by = join_by(Metacell == metacell)) %>%
                               rename(cell_type = cell.type) %>%
                               column_to_rownames(var = "Cell"))

# plot features and counts before cell filtering
mlei_beforeFiltering <- plot_featuresVScounts(mlei_S.object)

# plot features and counts before cell filtering
mlei_afterFiltering <- plot_featuresVScounts(mlei_S.object,
                                             nCount_min = 200, nCount_max = 4500,
                                             nFeature_min = 200, nFeature_max = 1300)

# convert SeuratObject to anndata file
fromSeurat_toAnndata(mlei_afterFiltering$S.object,
                     output_path, mlei_files$out_h5ad)


########################
#     N. VECTENSIS     #
########################

# create the Seurat object
nvec_S.object <- readRDS(nvec_files$RDS) %>%
  Seurat::CreateSeuratObject(project = "nvec", min.cells = 0, min.features = 0,
                             
                             # add meta data with pre-computed cell clusters and annotations
                             meta.data = read.table(nvec_files$metacell_assignments, header = TRUE, sep = "\t") %>%
                               left_join(read.table(nvec_files$metacell_annotation, header = TRUE, sep = "\t") %>%
                                           select(metacells, cell.type),
                                         by = join_by(Metacell == metacells)) %>%
                               rename(cell_type = cell.type) %>%
                               column_to_rownames(var = "Cell"))

# plot features and counts before cell filtering
nvec_beforeFiltering <- plot_featuresVScounts(nvec_S.object)

# plot features and counts before cell filtering
nvec_afterFiltering <- plot_featuresVScounts(nvec_S.object,
                                             nCount_min = 200, nCount_max = 5000,
                                             nFeature_min = 200, nFeature_max = 1500)

# convert SeuratObject to anndata file
fromSeurat_toAnndata(nvec_afterFiltering$S.object,
                     output_path, nvec_files$out_h5ad)


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

# convert SeuratObject to anndata file
fromSeurat_toAnndata(slac_afterFiltering$S.object,
                     output_path, slac_files$out_h5ad)


#########################
#     S. PISTILLATA     #
#########################

# create the Seurat object
spis_S.object <- readRDS(spis_files$RDS) %>%
  Seurat::CreateSeuratObject(project = "spis_adult", min.cells = 0, min.features = 0,
                             
                             # add meta data with pre-computed cell clusters and annotations
                             meta.data = read.table(spis_files$metacell_assignments, header = TRUE, sep = "\t") %>%
                               left_join(read.table(spis_files$metacell_annotation, header = TRUE, sep = "\t") %>%
                                           select(metacell, cell_type)) %>%
                               column_to_rownames(var = "cell"))

# plot features and counts before cell filtering
spis_beforeFiltering <- plot_featuresVScounts(spis_S.object)

# plot features and counts before cell filtering
spis_afterFiltering <- plot_featuresVScounts(spis_S.object,
                                             nCount_min = 200, nCount_max = 5000,
                                             nFeature_min = 200, nFeature_max = 1800)

# convert SeuratObject to anndata file
fromSeurat_toAnndata(spis_afterFiltering$S.object,
                     output_path, spis_files$out_h5ad)


########################
#     T. ADHAERENS     #
########################

# create the Seurat object fromt the UMI table
tadh_S.object <- load_placozoa_data(tadh_files$dir,
                                    tadh_files$mc,
                                    tadh_files$Rda,
                                    "tadh",
                                    tadh_files$metacell_annotation)

# plot features and counts before cell filtering
tadh_beforeFiltering <- plot_featuresVScounts(tadh_S.object)

# plot features and counts before cell filtering
tadh_afterFiltering <- plot_featuresVScounts(tadh_S.object,
                                             nCount_min = 200, nCount_max = 7500,
                                             nFeature_min = 200, nFeature_max = 2000)

# convert SeuratObject to anndata file
fromSeurat_toAnndata(tadh_afterFiltering$S.object,
                     output_path, tadh_files$out_h5ad)


###########################
#     T. ADHAERENS H2     #
###########################

# create the Seurat object fromt the UMI table
trh2_S.object <- load_placozoa_data(trh2_files$dir,
                                    trh2_files$mc,
                                    trh2_files$Rda,
                                    "trh2",
                                    trh2_files$metacell_annotation)

# plot features and counts before cell filtering
trh2_beforeFiltering <- plot_featuresVScounts(trh2_S.object)

# plot features and counts before cell filtering
trh2_afterFiltering <- plot_featuresVScounts(trh2_S.object,
                                             nCount_min = 200, nCount_max = 6000,
                                             nFeature_min = 200, nFeature_max = 1700)

# convert SeuratObject to anndata file
fromSeurat_toAnndata(trh2_afterFiltering$S.object,
                     output_path, trh2_files$out_h5ad)


##################
#     X. SP.     #
##################

# create the Seurat object
xesp_S.object <- readRDS(xesp_files$RDS) %>%
  Seurat::CreateSeuratObject(project = "xesp", min.cells = 0, min.features = 0,
                             
                             # add meta data with pre-computed cell clusters and annotations
                             meta.data = read.table(xesp_files$metacell_assignments, header = TRUE, sep = "\t") %>%
                               left_join(read.table(xesp_files$metacell_annotation, header = TRUE, sep = "\t") %>%
                                           select(metacell, cell_type)) %>%
                               column_to_rownames(var = "cell"))

# plot features and counts before cell filtering
xesp_beforeFiltering <- plot_featuresVScounts(xesp_S.object)

# plot features and counts before cell filtering
xesp_afterFiltering <- plot_featuresVScounts(xesp_S.object,
                                             nCount_min = 200, nCount_max = 7500,
                                             nFeature_min = 200, nFeature_max = 2800)

# convert SeuratObject to anndata file
fromSeurat_toAnndata(xesp_afterFiltering$S.object,
                     output_path, xesp_files$out_h5ad)
