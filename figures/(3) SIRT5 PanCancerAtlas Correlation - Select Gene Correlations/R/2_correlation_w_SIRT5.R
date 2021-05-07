#!/usr/bin/env Rscript

###---PACKAGES---###

if (!require("pacman")) { install.packages("pacman", repos='http://cran.us.r-project.org') }
library(pacman)

#required packages
required_packages = c(
  "tidyverse",
  "psych"
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
plot = "plot"

rsem.path = paste0(data.out,"SKCM_RNASeqV2.csv")

#load data
rsem = read.table(file = rsem.path, header = TRUE, sep = ',', check.names = FALSE) %>%
	column_to_rownames(var = "aliquot_barcode")

genes = names(rsem)
genes = genes[genes != "23408"] #exclude SIRT5 to avoid self correlation

cor.data.sirt5 = rsem %>%
	dplyr::select(c("23408")) %>%
	as.data.frame()
cor.data.genes = rsem %>%
	dplyr::select(genes) %>%
	as.data.frame()

spearman = corr.test(x = cor.data.sirt5, y = cor.data.genes, method = "spearman", adjust = "BH")
spearman.rho = t(spearman[[1]]) %>% as.data.frame() %>% setNames(c("spearman.rho")) %>% rownames_to_column(var="ENTREZID")
spearman.FDR = t(spearman[[4]]) %>% as.data.frame() %>% setNames(c("spearman.FDR")) %>% rownames_to_column(var="ENTREZID")

pearson = corr.test(x = cor.data.sirt5, y = cor.data.genes, method = "pearson", adjust = "BH")
pearson.cor = t(pearson[[1]]) %>% as.data.frame() %>% setNames(c("pearson.cor")) %>% rownames_to_column(var="ENTREZID")
pearson.FDR = t(pearson[[4]]) %>% as.data.frame() %>% setNames(c("pearson.FDR")) %>% rownames_to_column(var="ENTREZID")


res = spearman.rho %>% 
	left_join(spearman.FDR) %>% 
	left_join(pearson.cor) %>% 
	left_join(pearson.FDR)


write.csv(x = res, file = paste0(data.out,"correlation_w_SIRT5.csv"), row.names = FALSE)