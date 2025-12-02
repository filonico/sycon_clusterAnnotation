#!/usr/bin/env Rscript

setwd("/data/evassvis/fn76/sycon_clusterAnnotation")

library(Seurat)
library(tidyverse)
library(DESeq2)


#####################
#     LOAD DATA     #
#####################

load("./00_input/Sycon_Seuratv4.Rdata")

# get DE genes per cluster
# markers_deseq2 <- Sycon %>% FindAllMarkers(test.use = "DESeq2", verbose = TRUE, assay = "RNA", slot = "counts")
markers_wilcox <- Sycon %>% FindAllMarkers(test.use = "wilcox", verbose = TRUE)

# get the list of upregulated genes per cluster
markers_list <- markers_wilcox %>%
  filter(avg_log2FC > 0) %>%
  group_by(cluster) %>%
  # summarise(count = n())
  summarise(genes = list(gene), .groups = "drop") %>%
  { setNames(.$genes, .$cluster) }

# write theeach list of upregulated genes per cluster to a file
for (cluster_name in names(markers_list)) {
  file_path <- file.path("10_GO_enrichment", paste0("cluster", cluster_name, "_upregulatedGenes.ls"))
  writeLines(markers_list[[cluster_name]],
             file_path#,
             # row.names = FALSE,
             # col.names = FALSE,
             # quote = FALSE
            )
}

# write the list of gene universe to a file
writeLines(rownames(Sycon),
           file.path("10_GO_enrichment", "geneUniverse.ls"))
