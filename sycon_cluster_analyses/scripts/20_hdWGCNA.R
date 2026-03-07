#!/usr/bin/env Rscript

setwd("/data/evassvis/fn76/sycon/sycon_clusterAnnotation/sycon_cluster_analyses/")

suppressMessages({library(tidyverse)
  library(hdWGCNA)
  library(SeuratExtend)
  library(GOSemSim)
  library(smacof)})


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
               size = 1.5) +
    
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

plot_module_UMAP <- function(module_name, color_name) {
  
  umap_df <- Sycon@reductions$umap@cell.embeddings %>%
    as_tibble(rownames = NA) %>%
    rownames_to_column() %>%
    left_join(Sycon@meta.data %>% rownames_to_column()) %>%
    select(rowname, UMAP_1, UMAP_2, all_of(module_name)) %>%
    pivot_longer(-c(rowname, UMAP_1, UMAP_2),
                 names_to = "module",
                 values_to = "expr")
  
  # enforce symmetric range
  max_abs <- max(abs(umap_df$expr), na.rm = TRUE)
  
  umap_df <- umap_df %>%
    mutate(expr = pmax(pmin(expr, max_abs), -max_abs)) %>%
    arrange(expr)
  
  p <- ggplot(umap_df, aes(UMAP_1, UMAP_2)) +
    geom_point(aes(col = expr), size = 1.5) +
    
    scale_color_gradient2(low = "grey75", mid = "grey95", high = color_name,
                          midpoint = 0, limits = c(-max_abs, max_abs),
                          breaks = c(-max_abs, max_abs), labels = c("−", "+")) +
    
    ggtitle(module_name)
  
  return(p)
}


get_GOterm_semantic_mds <- function(goterm_enrich_file) {
  goerms_bp_elim <- read.table(goterm_enrich_file, sep = "\t", header = TRUE, quote = "") %>%
    filter(classicFisher < 0.05)
  
  annoDb_sycon <- read.table("09_gene_annotation/GOterms_OMA.tsv", sep = "\t", skip = 4) %>%
    select(V2, V5, V9) %>%
    mutate(V9 = case_when(V9 == "C" ~ "CC",
                          V9 == "P" ~ "BP",
                          V9 == "F" ~ "MF",
                          TRUE ~ "")) %>%
    rename(c("V2" = "GENE", "V5" = "GO", "V9" = "ONTOLOGY"))
  annoDb_sycon
  
  semantic_data <- godata(annoDb = annoDb_sycon, ont = "BP", computeIC = FALSE)
  
  semantic_similarity <- mgoSim(goerms_bp_elim$GO.ID, goerms_bp_elim$GO.ID,
                                semData = semantic_data, measure = "Wang", combine = NULL)
  
  semantic_dissimilarity <- sim2diss(semantic_similarity, method = "confusion", to.dist = TRUE)
  
  mds <- mds(semantic_dissimilarity)
  
  tibble_to_plot <- mds$conf %>%
    as_tibble(rownames = NA) %>%
    rownames_to_column(var = "GO") %>%
    left_join(goerms_bp_elim, by = join_by(GO == GO.ID)) %>%
    arrange(desc(classicFisher))
  
  return(tibble_to_plot)
}


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

md_arrows <- list(
  
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
  
  annotation_custom(grob = grid::textGrob(label = "MD1",
                                          x = unit(0, "mm"), y = unit(0, "mm") - unit(2.5, "mm"),
                                          just = c(0, 1), gp = grid::gpar(fontsize = 10))),
  
  annotation_custom(grob = grid::textGrob(label = "MD2",
                                          x = unit(0, "mm") - unit(2.5, "mm"), y = unit(0, "mm"),
                                          just = c(0, 0), rot = 90, gp = grid::gpar(fontsize = 10))),
  
  coord_cartesian(clip = "off"))


theme_for_plots <- theme(
  # aspect.ratio = 1,
  plot.background = element_rect(fill = "transparent", colour = NA), 
  panel.background = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_line(color = "grey90", lineend = "round"),
  panel.border = element_rect(colour = "black", linewidth = .6),
  legend.background = element_rect(fill = "transparent", colour = NA),
  legend.key = element_rect(fill = "transparent", colour = NA),
  legend.key.width = unit(.4, "cm"),
  legend.key.height = unit(.4, "cm"),
  legend.position = "right",
  legend.title = element_text(face = "bold"),
  # plot.title = element_text(size = 13, hjust = 0.0, vjust = 1.75, face = "bold"),
  axis.line = element_blank(),
  axis.ticks = element_line(colour = "black", linewidth = .4),
  axis.ticks.length = unit(0.10, "cm"),
  axis.text.x = element_text(color = "black",
                             margin = margin(t = 4, r = 0, b = 0, l = 0)),
  axis.text.y = element_text(color = "black",
                             margin = margin(t = 0, r = 4, b = 0, l = 0)),
  axis.title.y = element_text(angle = 90, size = 13,
                              margin = margin(t = 0, r = 10, b = 0, l = 0)),
  axis.title.x = element_text(angle = 0, size = 13,
                              margin = margin(t = 10, r = 0, b = 0, l = 0)),
  strip.text = element_text(color = "black", face = "bold", hjust = 0),
  strip.placement = "outside"
  # strip.background = element_rect(color = "black", linewidth = .6, linetype = "solid"),
  # strip.text.x = element_text(color = "black")
)


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

# compute eigengene-based connectivity (kME)
Sycon <- ModuleConnectivity(Sycon,
                            group.by = "seurat_clusters")


########################
#     PLOT MODULES     #
########################

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
  RNA_umap_list[[i]] <- plot_module_UMAP(i, i)
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


############################################
#     EXPORT LISTS OF GENE PER MODULES     #
############################################

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


#################################
#     UMAPS FOR PUBLICATION     #
#################################

SeuratExtend:::color_iwh("intense", n = 5)
#"#996742" "#5e824b" "#5d7880" "#7169a7" "#9f516c"
#"#c7956d" "#8fae6c" "#6db7ad" "#879dcf" "#cc88ad"
# "#b88340" "#9749ba" "#659d56" "#b64a56" "#757bad"

module1_umap <- plot_module_UMAP("yellow", "#a0520f") +
  labs(title = "Module signature") +
  
  theme_bw(base_size = 12) +
  theme_for_UMAPS +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")
module1_umap

module2_umap <- plot_module_UMAP("turquoise", "#af0fae") +
  theme_bw(base_size = 12) +
  theme_for_UMAPS +
  theme(plot.title = element_blank(),
        legend.position = "none")
module2_umap

module3_umap <- plot_module_UMAP("green", "#0faf19") +
  theme_bw(base_size = 12) +
  theme_for_UMAPS +
  theme(plot.title = element_blank(),
        legend.position = "none")
module3_umap

module4_umap <- plot_module_UMAP("brown", "#b90632") +
  theme_bw(base_size = 12) +
  theme_for_UMAPS +
  theme(plot.title = element_blank(),
        legend.position = "none")
module4_umap

module5_umap <- plot_module_UMAP("red", "#0f0fa0") +
  umap_arrows +
  theme_bw(base_size = 12) +
  theme_for_UMAPS +
  theme(plot.title = element_blank(),
        legend.position = "none",
        plot.margin = margin(0, 0, 6, 6, "mm"))
module5_umap


#######################################
#     ENRICHMENTS FOR PUBLICATION     #
#######################################

module1_semantic_mds <- get_GOterm_semantic_mds("12_hdWGCNA/01_RNA_assay/yellow_RNA_module_genes_GOterms_topGO_BP_elim.txt")
module2_semantic_mds <- get_GOterm_semantic_mds("12_hdWGCNA/01_RNA_assay/turquoise_RNA_module_genes_GOterms_topGO_BP_elim.txt")
module3_semantic_mds <- get_GOterm_semantic_mds("12_hdWGCNA/01_RNA_assay/green_RNA_module_genes_GOterms_topGO_BP_elim.txt")
module4_semantic_mds <- get_GOterm_semantic_mds("12_hdWGCNA/01_RNA_assay/brown_RNA_module_genes_GOterms_topGO_BP_elim.txt")
module5_semantic_mds <- get_GOterm_semantic_mds("12_hdWGCNA/01_RNA_assay/red_RNA_module_genes_GOterms_topGO_BP_elim.txt")

module1_goMD <-  module1_semantic_mds %>%
  ggplot(aes(D1, D2)) +
  geom_point(aes(size = Significant, colour = -log(classicFisher)), alpha = 0.7) +
  
  geom_label(data = . %>%
                     filter(GO %in% c("GO:0030833", "GO:0007015")) %>%
                     summarise(D1 = mean(D1), D2 = mean(D2), .groups = "drop"),
                   label = "Actin filament\norganization", fill = alpha("white", 0.8), size = 3,
                   label.padding = unit(0.3, "lines")) +

  geom_label(data = . %>%
                     filter(GO %in% c("GO:0030334", "GO:2000145", "GO:0040012", "GO:0040011")) %>%
                     summarise(D1 = mean(D1), D2 = mean(D2), .groups = "drop"),
                   label = "Cell\nmigration", fill = alpha("white", 0.8), size = 3,
                   label.padding = unit(0.3, "lines")) +

  geom_label(data = . %>%
                     filter(GO %in% c("GO:0007131", "GO:0035825", "GO:0140527", "GO:1903046", "GO:0140013", "GO:0061982", "GO:0007127", "GO:0051321")) %>%
                     summarise(D1 = mean(D1), D2 = mean(D2), .groups = "drop"),
                   label = "Meiotic\nprocesses", fill = alpha("white", 0.8), size = 3,
                   label.padding = unit(0.3, "lines")) +

  # geom_label_repel(data = . %>%
  #                    filter(GO %in% c("GO:0009168", "GO:0009127", "GO:0046125", "GO:0009120", "GO:0009157")) %>%
  #                    summarise(D1 = mean(D1), D2 = mean(D2), .groups = "drop"),
  #                  label = "Ribonucleoside\nmetabolism", fill = alpha("white", 0.8), size = 3,
  #                  label.padding = unit(0.3, "lines")) +
  
  # ggrepel::geom_text_repel(data = . %>%
  #                            filter(Significant >= 2),
  #                          aes(label = GO)) +
   
  labs(title = "Major GO semantics",
       x = "MD1", y = "MD2",
       colour = "-log(p-val)", size = "Sig. genes") +
  
  scale_size_continuous(range = c(1,14), breaks = c(5, 10)) +
  scale_color_gradient(high = "lightblue", low = "darkblue") +
  
  guides(size = guide_legend(override.aes = list(color = "#706db8")),
         color = guide_colorbar(barheight = 4, barwidth = 1)) +

  theme_bw(base_size = 12) +
  theme_for_UMAPS +
  theme(plot.margin = margin(6, 0, 0, 10, "mm"))

module1_goMD

module2_goMD <-  module2_semantic_mds %>%
  ggplot(aes(D1, D2)) +
  geom_point(aes(size = Significant, colour = -log(classicFisher)), alpha = 0.7) +
  
  geom_label(data = . %>%
               filter(GO %in% c("GO:0046942", "GO:0015711", "GO:0050877")) %>%
               summarise(D1 = mean(D1), D2 = mean(D2), .groups = "drop"),
             label = "Organic molecule\nmetabolism", fill = alpha("white", 0.8), size = 3,
             label.padding = unit(0.3, "lines")) +

  geom_label(data = . %>%
               filter(GO %in% c("GO:0032836","GO:0021554","GO:0097065","GO:0033333","GO:0033338","GO:0060039","GO:0048589","GO:0048468","GO:0048736",
                                "GO:0010720","GO:0048638","GO:0051094","GO:0072163","GO:0072164","GO:0001823","GO:0001657","GO:0060537","GO:0007399",
                                "GO:0010975","GO:0002064","GO:0072006","GO:0048731","GO:0007517","GO:0043588","GO:0050793", "GO:0048856")) %>%
               summarise(D1 = mean(D1), D2 = mean(D2), .groups = "drop"),
             label = "Anatomical\ndevelopment", fill = alpha("white", 0.8), size = 3,
             label.padding = unit(0.3, "lines")) +

  geom_label(data = . %>%
               filter(GO %in% c("GO:0007160", "GO0031589", "GO0042221", "GO0007155", "GO0030198", "GO0007156", "GO0007166")) %>%
               summarise(D1 = mean(D1), D2 = mean(D2), .groups = "drop"),
             label = "Cell adhesion\nand proliferation", fill = alpha("white", 0.8), size = 3,
             label.padding = unit(0.3, "lines")) +
  
  # ggrepel::geom_text_repel(data = . %>%
  #                            filter(Significant > 6),
  #                          aes(label = GO)) +
  
  labs(x = "MD1", y = "MD2",
       colour = "-log(p-val)", size = "Sig. genes") +
  
  scale_size_continuous(range = c(1,14), breaks = c(5, 10)) +
  scale_color_gradient(high = "lightblue", low = "darkblue") +
  
  guides(size = guide_legend(override.aes = list(color = "#706db8")),
         color = guide_colorbar(barheight = 4, barwidth = 1)) +
  
  theme_bw(base_size = 12) +
  theme_for_UMAPS +
  theme(plot.margin = margin(0, 0, 0, 10, "mm"))

module2_goMD

module3_goMD <-  module3_semantic_mds %>%
  ggplot(aes(D1, D2)) +
  geom_point(aes(size = Significant, colour = -log(classicFisher)),alpha = 0.7) +
  
  geom_label(data = . %>%
               filter(GO %in% c("GO:0007009", "GO:0071554", "GO:0006366")) %>%
               summarise(D1 = mean(D1), D2 = mean(D2), .groups = "drop"),
             label = "Cell membrane\norganization", fill = alpha("white", 0.8), size = 3,
             label.padding = unit(0.3, "lines")) +

  geom_label(data = . %>%
               filter(GO %in% c("GO:0007229", "GO:0007155")) %>%
               summarise(D1 = mean(D1), D2 = mean(D2), .groups = "drop"),
             label = "Cell adhesion", fill = alpha("white", 0.8), size = 3,
             label.padding = unit(0.3, "lines")) +

  geom_label(data = . %>%
               filter(GO %in% c("GO:0009913", "GO:0008544", "GO:0060429")) %>%
               summarise(D1 = mean(D1), D2 = mean(D2), .groups = "drop"),
             label = "Epithelium\ndevelopment", fill = alpha("white", 0.8), size = 3,
             label.padding = unit(0.3, "lines")) +
  
  # ggrepel::geom_text_repel(data = . %>%
  #                            filter(Significant >= 2),
  #                          aes(label = GO)) +
  
  labs(#title = "Enriched GOs",
    x = "MD1", y = "MD2",
    colour = "-log(p-val)", size = "Sig. genes") +
  
  scale_size_continuous(range = c(1,14), breaks = c(5, 10)) +
  scale_color_gradient(high = "lightblue", low = "darkblue") +
  
  guides(size = guide_legend(override.aes = list(color = "#706db8")),
         color = guide_colorbar(barheight = 4, barwidth = 1)) +
  
  coord_cartesian(clip = "off") +
  
  theme_bw(base_size = 12) +
  theme_for_UMAPS +
  theme(plot.margin = margin(0, 0, 0, 10, "mm"))

module3_goMD

module4_goMD <-  module4_semantic_mds %>%
  ggplot(aes(D1, D2)) +
  geom_point(aes(size = Significant, colour = -log(classicFisher)), alpha = 0.7) +

  geom_label(data = . %>%
               filter(GO %in% c("GO:0050801", "GO:0055080", "GO:0003169")) %>%
               summarise(D1 = mean(D1), D2 = mean(D2), .groups = "drop"),
             label = "Ion\nhomeostasis", fill = alpha("white", 0.8), size = 3,
             label.padding = unit(0.3, "lines")) +

  geom_label(data = . %>%
               filter(GO %in% c("GO:0034220", "GO:0008656", "GO:0002600", "GO0006820", "GO0005698")) %>%
               summarise(D1 = mean(D1), D2 = mean(D2), .groups = "drop"),
             label = "Ion transmembrane\ntransport", fill = alpha("white", 0.8), size = 3,
             label.padding = unit(0.3, "lines")) +

  geom_label(data = . %>%
               filter(GO %in% c("GO:0098609", "GO:0034446", "GO:0007160", "GO:0033627", "GO:0007156")) %>%
               summarise(D1 = mean(D1), D2 = mean(D2), .groups = "drop"),
             label = "Cell adhesion", fill = alpha("white", 0.8), size = 3,
             label.padding = unit(0.3, "lines")) +

  geom_label(data = . %>%
               filter(GO %in% c("GO:0009718", "GO:0071895")) %>%
               summarise(D1 = mean(D1), D2 = mean(D2), .groups = "drop"),
             label = "Mineralized tissue\ndevelopment", fill = alpha("white", 0.8), size = 3,
             label.padding = unit(0.3, "lines")) +

  geom_label(data = . %>%
               filter(GO %in% c("GO:0016055", "GO:0051493")) %>%
               summarise(D1 = mean(D1), D2 = mean(D2), .groups = "drop"),
             label = "Cytskeleton\norganization", fill = alpha("white", 0.8), size = 3,
             label.padding = unit(0.3, "lines")) +

  # ggrepel::geom_text_repel(data = . %>%
  #                            filter(Significant >= 2),
  #                          aes(label = GO)) +
  
  labs(#title = "Enriched GOs",
    x = "MD1", y = "MD2",
    colour = "-log(p-val)", size = "Sig. genes") +
  
  scale_size_continuous(range = c(1,14), breaks = c(5, 10)) +
  scale_color_gradient(high = "lightblue", low = "darkblue") +
  
  guides(size = guide_legend(override.aes = list(color = "#706db8")),
         color = guide_colorbar(barheight = 4, barwidth = 1)) +
  
  coord_cartesian(clip = "off") +
  
  theme_bw(base_size = 12) +
  theme_for_UMAPS +
  theme(plot.margin = margin(0, 0, 0, 10, "mm"))

module4_goMD

module5_goMD <-  module5_semantic_mds %>%
  ggplot(aes(D1, D2)) +
  geom_point(aes(size = Significant, colour = -log(classicFisher)), alpha = 0.7) +
  
  # geom_label(data = . %>%
  #              filter(GO %in% c("GO:0042989")) %>%
  #              summarise(D1 = mean(D1), D2 = mean(D2), .groups = "drop"),
  #            label = "Actin metabolism", fill = alpha("white", 0.8), size = 3,
  #            label.padding = unit(0.3, "lines")) +

  geom_label(data = . %>%
               filter(GO %in% c("GO:0006915", "GO:0007166", "GO:0043122", "GO:0071356")) %>%
               summarise(D1 = mean(D1), D2 = mean(D2), .groups = "drop"),
             label = "Cell adhesion\nand signalling", fill = alpha("white", 0.8), size = 3,
             label.padding = unit(0.3, "lines")) +
  
  # ggrepel::geom_text_repel(data = . %>%
  #                            filter(Significant >= 2),
  #                          aes(label = GO)) +
  
  labs(#title = "Enriched GOs",
    x = "MD1", y = "MD2",
    colour = "-log(p-val)", size = "Sig. genes") +
  
  scale_size_continuous(range = c(1,14), breaks = c(5, 10)) +
  scale_color_gradient(high = "lightblue", low = "darkblue") +
  
  guides(size = guide_legend(override.aes = list(color = "#706db8")),
         color = guide_colorbar(barheight = 4, barwidth = 1)) +
  
  md_arrows +
  
  theme_bw(base_size = 12) +
  theme_for_UMAPS +
  theme(plot.margin = margin(0, 0, 6, 10, "mm"))

module5_goMD


#################
#     PANEL     #
#################

panel_umaps <- ggpubr::ggarrange(module1_umap, module1_goMD,
                                 module2_umap, module2_goMD,
                                 module3_umap, module3_goMD,
                                 module4_umap, module4_goMD,
                                 module5_umap, module5_goMD,
                                 nrow = 5, ncol = 2, align = "h",
                                 widths = c(0.6,1),
                                 labels = "AUTO")
panel_umaps

ggsave("tmp.png",
       panel_umaps, device = "png",
       width = 12/1.5, height = 22/1.5, units = "in", dpi = 300, bg = "white")


# Sycon@reductions$umap@cell.embeddings %>%
#   as_tibble(rownames = NA) %>%
#   rownames_to_column() %>%
#   
#   # then add information from meta data, including the score for each module
#   left_join(Sycon@meta.data %>% rownames_to_column()) %>%
#   select(rowname, UMAP_1, UMAP_2, red, blue, yellow, green, turquoise, brown) %>%
#   
#   # pivot data longer
#   pivot_longer(-c(rowname, UMAP_1, UMAP_2), names_to = "module", values_to = "expr") %>%
#   arrange(expr) %>%
#   # group_by(module) %>%
#   
#   # plot data
#   ggplot(aes(UMAP_1, UMAP_2)) +
#   geom_point(aes(col = expr),
#              # alpha = 0.5, 
#              size = 0.7
#   ) +
#   
#   # define color gradient rules
#   # note that in the original hdWGCNA script they make the range symmetrical
#   # is it legit? Isn't this affecting visualization?
#   # here, we do not make the range symmetrical
#   # scale_color_gradient2(low = "grey75", mid = "grey95", high = module_name,
#   #                       midpoint = 0) +
#   labs(x = "UMAP 1", y = "UMAP 2",
#        color = "Module\nsignature", title = "Module 1") +
#   
#   annotation_custom(grob = grid::segmentsGrob(x0 = unit(0, "mm"), x1 = unit(12, "mm"),
#                                               y0 = unit(0, "mm"), y1 = unit(0, "mm"),
#                                               arrow = arrow(length = unit(2.5, "mm"),
#                                                             ends = "last", type = "open"),
#                                               gp = grid::gpar(col = "black", fill = "black", lwd = 1))) +
#   
#   annotation_custom(grob = grid::segmentsGrob(x0 = unit(0, "mm"), x1 = unit(0, "mm"),
#                                               y0 = unit(0, "mm"), y1 = unit(12, "mm"),
#                                               arrow = arrow(length = unit(2.5, "mm"),
#                                                             ends = "last", type = "open"),
#                                               gp = grid::gpar(col = "black", fill = "black", lwd = 1))) +
#   
#   annotation_custom(grob = grid::textGrob(label = "UMAP 1",
#                                           x = unit(0, "mm"), y = unit(0, "mm") - unit(2.5, "mm"),
#                                           just = c(0, 1), gp = grid::gpar(fontsize = 10))) +
#   
#   annotation_custom(grob = grid::textGrob(label = "UMAP 2",
#                                           x = unit(0, "mm") - unit(2.5, "mm"), y = unit(0, "mm"),
#                                           just = c(0, 0), rot = 90, gp = grid::gpar(fontsize = 10))) +
#   
#   coord_cartesian(clip = "off") +
#   
#   facet_wrap(. ~ module) +
#   
#   # make the plot squared
#   coord_fixed() +
#   
#   theme_bw(base_size = 12) +
#   theme_for_UMAPS +
#   theme(plot.title = element_text(hjust = 0),
#         plot.margin = margin(5, 5, 10, 8, "mm"))


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