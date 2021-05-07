#!/usr/bin/env Rscript

#used to produce files that are used by BART to generate output for 2_bart_plot.R

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

###---GLOBAL CONFIG---###
padj_threshold = 0.01
lfc_rna_threshold = 1

###---FUNCTIONS---###

data.in = "../data/in/dge/featureCounts_deseq2/table/result_lfcShrink/standardized/sirt5_kd_over_sirt5_nt/"
data.out = "../data/out/"

if(!dir.exists(data.out)) {
	dir.create(data.out)
}

#results
result_filtered = read.csv(paste0(data.in,"result-lfcShrink_stndrd-filt_anno-basic_padj1_lfc0.csv")) %>%
	dplyr::filter(external_gene_name != "SIRT5") %>%
	dplyr::filter(padj < padj_threshold) %>%
  dplyr::filter(abs(log2FoldChange) > lfc_rna_threshold)

#genes changed (excluding SIRT5, as we know it wasn't changed by a TF)
up_genes = result_filtered$external_gene_name[result_filtered$log2FoldChange > 0]
down_genes = result_filtered$external_gene_name[result_filtered$log2FoldChange < 0]

write.table(x=up_genes,file=paste0(data.out,"up_genes.txt"),row.names = FALSE, sep='\t', col.names = FALSE, quote=FALSE)
write.table(x=down_genes,file=paste0(data.out,"down_genes.txt"),row.names = FALSE,sep='\t',col.names = FALSE, quote=FALSE)


###---BART2---###

#Subsequently ran bart2-docker with the following commands:

#docker run --rm -ti -v "$PWD/data":/home/BARTv2.0/data/ -v "$PWD/bin":/home/BARTv2.0/bin/ -w /home/BARTv2.0/ bart2:latest /bin/bash
#bart2 geneset -i data/input/up_genes.txt -s hg38 --outdir data/output/
#bart2 geneset -i data/input/down_genes.txt -s hg38 --outdir data/output/