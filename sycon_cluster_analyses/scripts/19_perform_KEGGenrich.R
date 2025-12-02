#!/usr/bin/env Rscript

# setwd("/data/evassvis/fn76/sycon_clusterAnnotation")

suppressMessages({library(dplyr)
  library(clusterProfiler)})

args = commandArgs(trailingOnly = TRUE)

# args1 = gene universe with KEGG
# args2 = genes of interest (list)
# args3 = output name


######################
#     READ INPUT     #
######################

# read file with gene universe and KO annotation
universe_kegg_annotation <- read.table(args[1],
                                       sep = "\t", header = FALSE, fill = TRUE) %>%
  filter(V2 != "") %>%
  relocate(V2, .before = 1)

# read file of genes of interest
gene_interest <- read.table(args[2]) %>%
  unlist() %>%
  as.vector()


###################################
#     PERFORM KEGG ENRICHMENT     #
###################################

# get the list of KO terms for genes of interest
kegg_interest <- universe_kegg_annotation %>%
  filter(V1 %in% gene_interest) %>%
  select(V2) %>%
  unlist() %>%
  as.vector()

# get the list of KO terms for gene universe
kegg_universe <- universe_kegg_annotation %>%
  select(V2) %>%
  unlist() %>%
  as.vector()

# do KEGG enrichment
enrichKEGG <- enrichKEGG(kegg_interest,
                         organism = "ko",
                         keyType = "kegg",
                         universe = kegg_universe,
                         pvalueCutoff = 0.05)

# write results to file
write.table(x = enrichKEGG@result %>%
              filter(p.adjust <= 0.05) %>%
              arrange(category, p.adjust),
            file = args[3],
            quote = FALSE, row.names = FALSE, sep = "\t"
            )

