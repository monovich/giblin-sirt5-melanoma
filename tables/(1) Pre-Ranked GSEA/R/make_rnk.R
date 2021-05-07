#!/usr/bin/env Rscript

#RANK BASED GSEA
# 1) Generate adjusted Log2FC results (DESeq2)
# YOU ARE HERE --> 2) Generate rnk file for GSEA
# 3) GSEA tool analysis

###---PACKAGES---###

if (!require("pacman")) { install.packages("pacman", repos='http://cran.us.r-project.org') }
library(pacman)

#required packages
required_packages = c(
  "tidyverse"
)

#load packages
pacman::p_load(
  char=required_packages,
  install=TRUE,
  character.only=TRUE,
  try.bioconductor=TRUE,
  update.bioconductor=TRUE
)

###---RNK---###

data.in = "../data/in/dge/featureCounts_deseq2/table/result/standardized/sirt5_kd_over_sirt5_nt/"
data.out = "../data/out/"

src = paste0(data.in,"result-lfcShrink_stndrd-filt_anno-basic_padj1_lfc0.csv")

rnk.src = read.csv(src) %>%
  as_tibble()

rnk.sort = rnk.src %>%
  dplyr::select(c("external_gene_name","log2FoldChange")) %>%
  dplyr::arrange(-log2FoldChange)

write.table(
  x = rnk.sort,
  file = paste0(data.out,"sirt5_kd_over_sirt5_nt.rnk"),
  row.names = F,
  col.names = F,
  quote = F,
  sep = "\t"
)

###---GSEAPreRanked---###

#INSTRUCTIONS:

#1) Load GSEA Application (tested with GSEA 4.1.0 ; MsigDB v7.2)
#2) Load .rnk file via "Load Data" in lefthand menu
#3) Run "GSEAPreRanked" Analysis on loaded .rnk file with following parameters:
# Required Fields
# - Gene sets database (one at a time):
#   - h.all.v7.2.symbols.gmt [Hallmarks]
#   - c1.all.v7.2.symbols.gmt [Positional]
#   - c2.all.v7.2.symbols.gmt [Curated]
#   - c3.all.v7.2.symbols.gmt [Motif]
#   - c4.all.v7.2.symbols.gmt [Computational]
#   - c5.go.bp.v7.2.symbols.gmt [Gene ontology]
#   - c5.go.cc.v7.2.symbols.gmt [Gene ontology]
#   - c5.go.mf.v7.2.symbols.gmt [Gene ontology]
#   - c6.all.v7.2.symbols.gmt [Oncogenic signatures]
#   - c7.all.v7.2.symbols.gmt [Immunologic signatures]
#   - c8.all.v7.2.symbols.gmt [Cell type signatures]

# - Number of permutations: 1000

# - Ranked list: loaded .rnk file

# - Collapse/Remap to gene symbols: Collapse

# - Chip platform: Human_Gene_Symbol_with_Remapping_MSigDB.v7.2.chip

# Basic Fields (default options)
# - Enrichment statistic: weighted

# - Max size: 500

# - Min size: 15

# Advanced fileds
# - Collapsing mode for probe sets => 1 gene: Max_probe
# - Normalization mode: meandiv
# - Create SVG plot images: true
# - Omit features with no symbol match: true
# - Make detailed gene set report: true
# - Plot graphs for the top sets of each phenotype: 20
# - Seed permutation: timestamp
# - Make a zipped file with all reports: true

