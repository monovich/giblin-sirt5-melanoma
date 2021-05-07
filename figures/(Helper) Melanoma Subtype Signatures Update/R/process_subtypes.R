#!/usr/bin/env Rscript

#NOTES:

#tsoi_2018_supplement data taken directly from:
# Tsoi, J., Robert, L., Paraiso, K., Galvan, C., Sheu, K. M., Lay, J., ... & Graeber, T. G. (2018). 
#   Multi-stage differentiation defines melanoma subtypes with differential vulnerability to drug-induced 
#   iron-dependent oxidative stress. Cancer cell, 33(5), 890-904.

###---PACKAGES---###

if (!require("pacman")) { install.packages("pacman", repos='http://cran.us.r-project.org') }
library(pacman)

#required packages
required_packages = c(
  "tidyverse",
  "readxl"
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

data.input = "data/input/"
data.output = "data/output/"

subtype_signatures_original = read_xlsx(
	skip=1,
	path = paste0(data.input,"tsoi_2018_supplement/1-s2.0-S1535610818301223-mmc4.xlsx")
)

write.csv(
	x = subtype_signatures_original %>%
		dplyr::select(c("Gene","Signature")),
	file = paste0(data.output,"subtype_signatures_original.csv"),
	row.names = FALSE
)

### RUN FILE THROUGH HUGO MULTI-SYMBOL CHECKER
# https://www.genenames.org/tools/multi-symbol-checker/
# for input subtype_signatures_original.csv, output file manually named subtype_signatures_updated.csv

subtype_signatures_updated = read.csv(file = paste0(data.output,"subtype_signatures_updated.csv"))