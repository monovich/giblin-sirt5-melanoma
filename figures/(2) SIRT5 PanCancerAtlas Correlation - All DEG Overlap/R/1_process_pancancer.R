#!/usr/bin/env Rscript

#PanCancerAtlas data retrieved from: https://gdc.cancer.gov/about-data/publications/pancanatlas

###---PACKAGES---###

if (!require("pacman")) { install.packages("pacman", repos='http://cran.us.r-project.org') }
library(pacman)

#required packages
required_packages = c(
  "tidyverse",
  "data.table"
)

#load packages
pacman::p_load(
  char=required_packages,
  install=TRUE,
  character.only=TRUE,
  try.bioconductor=TRUE,
  update.bioconductor=TRUE
)

###---FUNCTIONS---###

data.in = "data/in/"
data.out = "data/out/"

sample_lookup = fread(
	file = paste0(data.in,"merged_sample_quality_annotations.tsv"),
	data.table = FALSE
)

sample_ids = sample_lookup$aliquot_barcode[sample_lookup[["cancer type"]] == "SKCM" & sample_lookup[["platform"]] == "IlluminaHiSeq_RNASeqV2"]

pancan = fread(
	file=paste0(data.in,"EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv"),
	select = c("gene_id",sample_ids),
	data.table = FALSE
)

pos1 = function(x) {
	return(unlist(strsplit(x = x, split = '|', fixed=TRUE))[1])
}

pos2 = function(x) {
	return(unlist(strsplit(x = x, split = '|', fixed=TRUE))[2])
}

pancan$entrez_id = sapply(pancan$gene_id, FUN = pos2)
pancan$gene_name = sapply(pancan$gene_id, FUN = pos1)

#SIRT5 (23408)

matrix = pancan %>% dplyr::select(c(sample_ids,"entrez_id")) %>%
	column_to_rownames(var="entrez_id")

matrix = t(matrix) %>%
	as.data.frame() %>%
	rownames_to_column(var="aliquot_barcode")

write.csv(x = matrix, file = paste0(data.out, "SKCM_RNASeqV2.csv"),row.names = FALSE)