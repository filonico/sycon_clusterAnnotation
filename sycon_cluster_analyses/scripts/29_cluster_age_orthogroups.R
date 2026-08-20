#!/usr/bin/env Rscript

setwd("/data/evassvis/fn76/sycon/sycon_clusterAnnotation/sycon_cluster_analyses/")

library(tidyverse)
library(Seurat)
library(SeuratExtend)
library(patchwork)


###########################################
#     LOAD ORTHOGROUP NODE ASSIGNMENT     #
###########################################

# create tibble to store og-to-node mapping
og_node_assignement <- tibble("og" = character(),
                              "value" = double(),
                              "node" = character())

# populate tibble with node files
for (node in seq(116, 230)) {
  # read in the file
  node_file <- read.table(paste0("18_cluster_age_orthogroups/01_node_lists/", node, ".txt"),
                          header = FALSE, col.names = c("og", "value"), sep = "\t") %>%
    # add node name column
    mutate(node = paste0("N", node))
  # append to tibble
  og_node_assignement <- og_node_assignement %>% add_row(node_file)
  
  rm(node_file)
}
og_node_assignement

# read in table with nodes of interest
nodes_of_interest <- read.table("18_cluster_age_orthogroups/notable_nodes.tsv",
                                header = FALSE, col.names = c("node", "clade"),
                                sep = "\t")
nodes_of_interest
