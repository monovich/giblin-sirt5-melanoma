#!/usr/bin/env Rscript

# DESeq2 analysis script

############################################
#            PACKAGES
############################################

#high usage
if (!require("pacman")) { install.packages("pacman") }
library(pacman)

#required packages
required_packages = c(
  "tidyverse",
  "DESeq2",
  #"fgsea",
  #"msigdbr",
  "BiocParallel",
  "GenomicFeatures",
  "ensembldb",
  "optparse",
  "apeglm"
)

github_packages = c(
  "jmw86069/splicejam"
)

#load packages
pacman::p_load(
  char = required_packages,
  install = TRUE,
  force = FALSE,
  character.only = TRUE,
  try.bioconductor = TRUE,
  update.bioconductor = TRUE
)

#load github packages
pacman::p_load_gh(
  char = github_packages,
  update = getOption("pac_update"),
  dependencies = TRUE
)

# ############################################
# #            ARGUMENTS
# ############################################

#sample execution:
# Rscript -p ~/project -d ~/project/data -g ~/genome -e mmusculus_gene_ensembl

# Optparser
option_list = list(
  make_option(c("-p", "--projectDir"),
              action="store",
              type="character",
              default=getwd(), 
              help="working directory path",
              metavar="character"),
  make_option(c("-d", "--dataDir"),
              action="store",
              type="character",
              default=NULL, 
              help="data directory path",
              metavar="character"),
  make_option(c("-g", "--genomeDir"),
              action="store",
              type="character",
              default=NULL, 
              help="genome directory path",
              metavar="numeric"),
  make_option(c("-e", "--ensemblSpecies"),
              action="store",
              type="character",
              default="hsapiens_gene_ensembl", 
              help="ensembl species dataset (e.g. hsapiens_gene_ensembl, mmusculus_gene_ensembl)",
              metavar="character"),
  make_option(c("-t", "--threads"),
              action="store",
              type="numeric",
              default=parallel::detectCores()-2, 
              help="threads available",
              metavar="numeric")
  # make_option(c("-s", "--contrastsToStandardizeBy"),
  #             action="store",
  #             type="character",
  #             default=NULL, 
  #             metavar="character")
) 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

genome_dir = paste0(opt$g,"/")
workers_available = opt$t
ensembl_dataset = opt$e
#standardizedContrasts= as.numeric(unlist(strsplit(split=",",x=opt$s)))
options(max.print=100) #for sanity

############################################
#            TESTING
############################################

path="/run/media/monovich/files/phd/billy_sirt5_kd/"
project=paste0(path,"final_versions/DESeq2/")
genome=paste0(path,"genome/")
data=paste0(project,"data/")

opt=list()
genome_dir = paste0(path,"/genome/gencode/Gencode_human/GRCh38/release_34/")
opt$p = project
opt$dataDir = data
workers_available = parallel::detectCores()-2
ensembl_dataset = "hsapiens_gene_ensembl"
options(max.print=100)

############################################
#      Initialize Master Directories/Files
############################################

#important master directories
data_in = paste0(opt$dataDir,"/in/")
data_out = paste0(opt$dataDir,"/out/")

#important sub directories
featureCounts_out = paste0(data_out,"featureCounts/trimmed/paired/")
dge_in = paste0(data_in,"/dge/")
dge_out = paste0(data_out,"/dge/")
featureCounts_deseq2_out = paste0(dge_out,"/featureCounts_deseq2/")

#important files
annotation_gtf = paste0(genome_dir,"annotation.gtf")
txdb = makeTxDbFromGFF(file=annotation_gtf)
exons_by_gene = exonsBy(txdb, by = "gene")

############################################
#       Fetch Feature Counts Matrices
############################################

#featureCounts_out="/run/media/monovich/files/phd/billy_sirt5_kd/project/sirt5kd_mel_indiv_allsamples_unmerged/data/out/featureCounts/trimmed/paired/"

#fetch all featureCounts txt matrix paths
featureCounts_matrices = list.files(
  path=featureCounts_out,
  full.names = TRUE,
  pattern="*.txt$"
)

#get sample ids
featureCounts_sample_ids = basename(gsub(pattern=".txt","",featureCounts_matrices))

#add names to paths
names(featureCounts_matrices) = featureCounts_sample_ids

#usage: featureCounts_matrices['87752']

############################################
#       Generate Feature Annotations
############################################

annotation_out = paste0(genome_dir,"annotation/")
if (!dir.exists(annotation_out)) { dir.create(annotation_out,recursive = T) }

gene_feature_annotations.csv = paste0(annotation_out, "gene_feature_annotations.csv")

if (! file.exists(gene_feature_annotations.csv)) {
  
  Tx2geneFromGtf.csv = paste0(annotation_out,"Tx2geneFromGtf.csv")
  
  if (! file.exists(Tx2geneFromGtf.csv)) {
            
    #make transcript_feature_type for identifying gene_feature_type
    Tx2geneFromGtf = splicejam::makeTx2geneFromGtf(annotation_gtf)
    write.csv(x=Tx2geneFromGtf,
              file=gene_feature_annotations.csv,
              row.names = F)
  } else {
    Tx2geneFromGtf = read.csv(file=Tx2geneFromGtf.csv)
  }

  #filter to gene level
  gene_feature_annotations = Tx2geneFromGtf %>%
    dplyr::select(c("gene_id","gene_type")) %>% 
    distinct() %>% 
    setNames(c("ensembl_gene_id_version","gene_type"))
  
  #write annotation
  write.csv(x=gene_feature_annotations,
            file=gene_feature_annotations.csv,
            row.names = F)
} else {
  #load annotations
  gene_feature_annotations = read.csv(file=gene_feature_annotations.csv,header=TRUE)
}

############################################
#            Import & Clean sampleLookup
############################################

sampleLookup = read.csv(file=paste0(dge_in,"sampleLookup.csv"),header=TRUE) 

############################################
#            Import groupContrast
############################################

#read Contrasts
groupContrast = read.csv(file=paste0(dge_in,"groupContrast.csv"),header=TRUE)

#all models/contrasts specified (no all/~1)
contrasts = groupContrast$contrast
models = groupContrast$model

#create list of sample ids to subset by for comparisons, named by contrast
contrastSampleIds = list()
contrastTxtPaths = list()
contrastColData = list()

model_to_formula = function(x) {
  return(as.formula(paste0("~",x)))
}

#ADD CONTRASTS INFO
for (i in 1:length(groupContrast$contrast)) {
  top_colData = sampleLookup
  top_cols = unlist(strsplit(x = groupContrast$top_cols[i],split = "+",fixed = TRUE))
  top_vals = unlist(strsplit(x = groupContrast$top_vals[i],split = "+",fixed = TRUE))
  
  for(j in 1:length(top_cols)) {
    top_colData = subset(top_colData, top_colData[[top_cols[j]]] %in% top_vals[j])
  }
  #add unique_group_id
  top_colData = top_colData %>% 
    mutate(unique_group_id=as.factor(groupContrast$top_vals[i]))
  
  bottom_colData = sampleLookup
  bottom_cols = unlist(strsplit(x = groupContrast$bottom_cols[i],split = "+",fixed = TRUE))
  bottom_vals = unlist(strsplit(x = groupContrast$bottom_vals[i],split = "+",fixed = TRUE))
  
  for(j in 1:length(bottom_cols)) {
    bottom_colData = subset(bottom_colData, bottom_colData[[bottom_cols[j]]] %in% bottom_vals[j])
  }
  #add unique_group_id
  bottom_colData = bottom_colData %>% 
    mutate(unique_group_id=as.factor(groupContrast$bottom_vals[i]))
  
  #model terms
  model.terms = unique(unlist(strsplit(x = str_remove(groupContrast$model[i],"\\~"),split = "[+:]",fixed = F)))
  
  #levels
  contrast_levels = list()
  tmp = unlist(strsplit(x = groupContrast$levels[i],split = "[+]",fixed = F))
  for (j in 1:length(tmp)) {
    contrast_levels[[j]] = unlist(strsplit(x = tmp[j],split = "[,]",fixed = F))
  }
  
  #arrange terms in reverse order (ensures last var of interest is dominant sort)
  arrange = rev(unlist(strsplit(x = groupContrast$arrange[i],split = "+",fixed = T)))
  
  #rbind colData, control rows first, treatment second
  contrastColData[[i]] = rbind(bottom_colData,top_colData) %>% 
    arrange_(.dots=arrange)
  
  #factor terms by predetermined order
  for(j in 1:length(model.terms)) {
    print("Before Factoring:")
    print(contrastColData[[i]][[model.terms[j]]])
    contrastColData[[i]][[model.terms[j]]] = factor(x=contrastColData[[i]][[model.terms[j]]],
                                                    levels=contrast_levels[[j]])
    print("After Factoring:")
    print(contrastColData[[i]][[model.terms[j]]])
  }
  
  #provide row names that just equal sample IDs
  row.names(contrastColData[[i]]) = contrastColData[[i]]$sample_id
  
  #store sample IDs and paths
  contrastSampleIds[[i]] = c(bottom_colData$sample_id,top_colData$sample_id)
  contrastTxtPaths[[i]] = subset(featureCounts_matrices, names(featureCounts_matrices) %in% contrastSampleIds[[i]])
  
  #makes matrix column order match sample_id row order as columns are added by file order
  contrastTxtPaths[[i]] = contrastTxtPaths[[i]][order(factor(names(contrastTxtPaths[[i]]), levels = contrastColData[[i]]$sample_id))]
}

#ADD FULL TABLE INFO (~1)
contrastSampleIds[[1+length(contrastSampleIds)]] = sampleLookup$sample_id
contrastTxtPaths[[1+length(contrastTxtPaths)]] = subset(featureCounts_matrices, names(featureCounts_matrices) %in% sampleLookup$sample_id)
contrastColData[[1+length(contrastColData)]] = sampleLookup

#UPDATE MODELS
models=c(models,"1")
contrasts=c(contrasts,"all")

#ADD NAMES
names(contrastSampleIds) = contrasts
names(contrastTxtPaths) = contrasts
names(contrastColData) = contrasts


############################################
#            Make Model Matrices
############################################

contrastModelMatrix = list()
for (i in 1:length(models)) {
  contrastModelMatrix[[i]] = model.matrix(object=model_to_formula(models[i]),data=contrastColData[[i]])
}
# contrastModelMatrix[[length(contrastModelMatrix)+1]] = model.matrix(object=~1,data=sampleLookup)
names(contrastModelMatrix) = contrasts

############################################
#       Initialize Plot Directories
############################################

#plot output directories
plot.dir = paste0(featureCounts_deseq2_out,"/plot/")
plot.dirs = list()
for (i in 1:length(contrasts)) {
  contrast = contrasts[i]
  contrast = iconv(contrast, from = 'UTF-8', to = 'ASCII//TRANSLIT')
  contrast = gsub("?","a",contrast,fixed=T)
  plot.dirs[[i]] = paste0(plot.dir,"/",contrast,"/")
}
plot.dirs = unlist(plot.dirs)

############################################
#       Read FeatureCounts Matrices
############################################

count_matrices = list()
for (i in 1:length(contrasts)) {
  ordered_sample_ids = contrastColData[[i]]$sample_id

  for (j in 1:length(ordered_sample_ids)) {
    if (j == 1) {
      id = toString(ordered_sample_ids[j])
      tmp_matrix = read.csv(featureCounts_matrices[id],
                             sep="",
                             head=TRUE,
                             skip=1,
                             row.names = "Geneid")
      count_matrix = as.data.frame(tmp_matrix[,6])
      row.names(count_matrix) = row.names(tmp_matrix)
      colnames(count_matrix) = id
      count_matrix = count_matrix %>% 
        rownames_to_column()
      count.matrix = count_matrix
    } else {
      id = toString(ordered_sample_ids[j])
      tmp_matrix = read.csv(featureCounts_matrices[id],
                            sep="",
                            head=TRUE,
                            skip=1,
                            row.names = "Geneid")
      count_matrix = as.data.frame(tmp_matrix[,6])
      row.names(count_matrix) = row.names(tmp_matrix)
      colnames(count_matrix) = id
      count_matrix = count_matrix %>% 
        rownames_to_column()
      count.matrix = count.matrix %>% 
        full_join(count_matrix, by="rowname")
    }
  }

  count.matrix = count.matrix %>% column_to_rownames("rowname")
  count_matrices[[i]] = count.matrix
}
names(count_matrices) = contrasts

############################################
#       DDS Construction
############################################

#Configure multicore usage
register(MulticoreParam(workers=workers_available, exportglobals = FALSE))

#function for getting ddsMatrix
get_dds=function(countData,colData,contrast,formula.string,count.threshold=10,sample.threshold=3,mean.threshold=1) {
  print(paste0(contrast))
  ddsMatrix = DESeqDataSetFromMatrix(countData = countData,
                                     colData = colData,
                                     design = model_to_formula(formula.string)) 
  #process
  if(contrast == "all") {
    ddsMatrix = DESeq(ddsMatrix, parallel = TRUE, BPPARAM=MulticoreParam(workers=workers_available,progressbar=F),minReplicatesForReplace=Inf)
  } else {
    ddsMatrix = DESeq(ddsMatrix, parallel = TRUE, BPPARAM=MulticoreParam(workers=workers_available,progressbar=F))
  }
  
  #Remove rows that don't converge in beta
  ddsMatrix[which(mcols(ddsMatrix)$betaConv),]
  
  return(ddsMatrix)
}

ddsMatrices = list()
for (i in 1:(length(contrasts)-1)) {
  ddsMatrices[[i]]= get_dds(
    countData=count_matrices[[i]],
    colData=contrastColData[[i]],
    contrast=contrasts[i],
    formula.string=models[i]
  )
}

ddsMatrix_all = get_dds(
  countData=count_matrices[["all"]],
  colData=contrastColData[["all"]],
  contrast="all",
  formula.string="1"
)

############################################
#      Fetch Annotations
############################################

#make dir
ensembl.dir = paste0(annotation_out,"/ensembl/")
if (!dir.exists(ensembl.dir)) { dir.create(ensembl.dir,recursive = T) }

#gene annotations
if (! file.exists(paste0(ensembl.dir,"gene_annotations.csv"))) {
  mart = biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL",
                          dataset=ensembl_dataset, 
                          host = 'uswest.ensembl.org')
  attributes=c(
    "external_gene_name",
    "description",
    "ensembl_gene_id",
    "ensembl_gene_id_version"
  )
  gene_annotations = biomaRt::getBM(attributes=attributes,mart=mart)
  write.csv(gene_annotations,paste0(ensembl.dir,"gene_annotations.csv"),row.names = F)
}

gene_annotations = read.csv(file=paste0(ensembl.dir,"gene_annotations.csv"))

#gene annotations to homolog
if (ensembl_dataset != "hsapiens_gene_ensembl") {
  if (! file.exists(paste0(ensembl.dir,"homolog_annotations.csv"))) {
    mart = biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL",
                            dataset=ensembl_dataset, 
                            host = 'uswest.ensembl.org')
    attributes=c(
      "hsapiens_homolog_associated_gene_name",
      "ensembl_gene_id_version"
    )
    homolog_annotations = biomaRt::getBM(attributes=attributes,mart=mart) %>% 
      setNames(c("human_gene_symbol","ensembl_gene_id_version")) %>%
      distinct() %>%
      as_tibble() %>%
      na_if("") %>% 
      na.omit()
    write.csv(homolog_annotations,paste0(ensembl.dir,"homolog_annotations.csv"),row.names = F)
  }

  homolog_annotations = read.csv(file=paste0(ensembl.dir,"homolog_annotations.csv"))
}

# lists of annotations to merge
ensembl = list(
  gene_annotations,
  gene_feature_annotations
)

if (ensembl_dataset != "hsapiens_gene_ensembl") {
  ensembl_homolog = list(
    gene_annotations,
    homolog_annotations,
    gene_feature_annotations
  )
}

if (ensembl_dataset != "hsapiens_gene_ensembl") {
  ensembl_entrez_homolog = list(
    gene_annotations,
    homolog_annotations,
    gene_feature_annotations
  )
}

if (ensembl_dataset != "hsapiens_gene_ensembl") {
  possible_annotations = list(
    ensembl,
    ensembl_homolog
  )
} else {
  possible_annotations = list(
    ensembl
  )
}

if (ensembl_dataset != "hsapiens_gene_ensembl") {
  names(possible_annotations) = c(
    "anno-basic",
    "anno-homolog"
  )
} else {
  names(possible_annotations) = c(
    "anno-basic"
  )
}

############################################
#    Get Filtered Results + Filtered DDS
############################################

#make dir
if (!dir.exists(featureCounts_deseq2_out)) { dir.create(featureCounts_deseq2_out,recursive = T) }
table_dir = paste0(featureCounts_deseq2_out,"/table/")
if (!dir.exists(table_dir)) { dir.create(table_dir,recursive = T) }

#getResults
results_independent = list()
for (i in 1:(length(ddsMatrices))) {
  print(paste0("results : independent (v1) : ",contrasts[i]))
  results_independent[[i]] = results(
    object = ddsMatrices[[i]],
    parallel = TRUE, 
    BPPARAM = SnowParam(workers = workers_available, progressbar = F),
    alpha = 0.05,
    independentFiltering = TRUE
  )
}

#getResults no indpendent filtering
results_nofilter = list()
for (i in 1:(length(ddsMatrices))) {
  print(paste0("results : no filter : ",contrasts[i]))
  results_nofilter[[i]] = results(
    object = ddsMatrices[[i]],
    parallel = TRUE, 
    BPPARAM = SnowParam(workers = workers_available, progressbar = F),
    independentFiltering = FALSE
  )
}

#generate a shared list of genenames that are not independently filtered
genenames = sapply(
  X = results_independent, 
  FUN = function(x) { rownames(as.data.frame(x))[!is.na(as.data.frame(x)$padj)] }
)

# genenames = list()
# getSharedGenenames = 
# j = 1
# for (i in standardizedContrasts) {
#   genenames[[j]] = getSharedGenenames(results_independent[[i]])
#   j = 1 + j
# }

combined_filtered_genenames = do.call(c, genenames) %>% unique()
write.table(x=combined_filtered_genenames,file=paste0(table_dir,"standardized_rownames.csv"),sep=",",row.names=FALSE,col.names=FALSE)

#apply filts to dds independently
ddsMatrices_independent=list()
for(i in 1:(length(ddsMatrices))) {
  keep = row.names(counts(ddsMatrices[[i]])) %in% genenames[[i]]
  ddsMatrices_independent[[i]] = ddsMatrices[[i]][keep,]
}

#standardize matrices filters
ddsMatrices_standardized=list()
for (i in 1:(length(ddsMatrices))) {
  keep = row.names(counts(ddsMatrices[[i]])) %in% combined_filtered_genenames
  ddsMatrices_standardized[[i]] = ddsMatrices[[i]][keep,]
}

#all
keep = row.names(counts(ddsMatrices[[i]])) %in% combined_filtered_genenames
ddsMatrix_all_standardized = ddsMatrix_all[keep,]

#regenerate results to get consistent rownames
#getResults
results_independent = list()
for (i in 1:(length(ddsMatrices_independent))) {
  print(paste0("results : independent (v2) : ",contrasts[i]))
  results_independent[[i]] = results(
    object = ddsMatrices_independent[[i]],
    parallel = TRUE, 
    BPPARAM = SnowParam(workers = workers_available, progressbar = F),
    independentFiltering = FALSE
  )
}

#export results with standardized filters
results_standardized = list()
for (i in 1:(length(ddsMatrices_standardized))) {
  print(paste0("results : standardized : ",contrasts[i]))
  results_standardized[[i]] = results(
    object = ddsMatrices_standardized[[i]],
    parallel = TRUE, 
    BPPARAM = SnowParam(workers = workers_available, progressbar = F),
    independentFiltering = FALSE
  )
}

############################################
#      Exports Preconfig
############################################

### export options
blind_options = c(FALSE) #vst,rlog
robust_options = c(TRUE) #fpkms
normalized_options = c(TRUE) #counts
padj_thresholds = c(1,0.05,0.01)
abs_lfc_thresholds = c(0,0.58)
methods = c("rlog","fpkm","count")
filters = c("independent","standardized")


ddsMatrices_complete=list(ddsMatrices,ddsMatrices_independent,ddsMatrices_standardized)
names(ddsMatrices_complete) = filters

ddsMatrices_all_complete=list(ddsMatrix_all,ddsMatrix_all_standardized)
names(ddsMatrices_all_complete) = c("none","standardized")

results_complete=list(results_nofilter,results_independent,results_standardized)
names(results_complete) = filters

############################################
#      Export Results
############################################

#function
exportResults = function(
  result,
  contrast,
  annotations,
  annotations.names,
  table.dir,
  filter=c("none","independent","standardized"),
  padj.thresholds=c(1),
  abs.lfc.thresholds=c(0)
) {

  #join by
  col = "ensembl_gene_id_version"

  #filter_suffix
  filter = match.arg(filter)
  filter_suffix = switch(filter,
    none="no-filt",
    independent="indep-filt",
    standardized="stndrd-filt"
  )

  #result_df
  result_df = as.data.frame(result) %>%
    rownames_to_column(col)

  #anno
  results_df_anno=list()
  for (i in 1:length(annotations)) {
    if (!is.null(annotations[[i]])) {
      tmp = result_df
      for (j in 1:length(annotations[[i]])) {
        tmp = tmp %>% left_join(annotations[[i]][[j]], by=col)
      }
      results_df_anno[[i]] = tmp
    } else {
      results_df_anno[[i]] = result_df
    }
  }

  #filter result_df for only significant genes
  result_df_signif = list()
  filename = list()

  l = 1
  for(i in 1:length(results_df_anno)) {
    for (j in 1:length(padj.thresholds)) {
      for (k in 1:length(abs.lfc.thresholds)) {
        result_df_signif[[l]] = results_df_anno[[i]] %>%
          dplyr::filter(padj < padj.thresholds[j]) %>%
          dplyr::filter(abs(log2FoldChange) > abs.lfc.thresholds[k])
        signif_suffix = paste0("padj",padj.thresholds[j],"_lfc",abs.lfc.thresholds[k])
        filename[[l]] = paste0("result","_",filter_suffix,"_",annotations.names[i],"_",signif_suffix,".csv")
        l = l + 1
      }
    }
  }
  filename=unlist(filename)

  #make dir
  final.dir = paste0(table.dir,"/",contrast,"/")
  if(!dir.exists(final.dir)) { dir.create(final.dir,recursive = T) }

  #write tables
  for (i in 1:length(result_df_signif)) {
    write.csv(x=result_df_signif[[i]], file=paste0(final.dir,filename[i]), row.names = F)
  }
  print(contrast)
  print(filename)
}

for (i in 1:(length(contrasts)-1)) {
  for (j in filters) {
    exportResults(
      result=results_complete[[j]][[i]],
      contrast=contrasts[i],
      annotations=possible_annotations,
      annotations.names=names(possible_annotations),
      table.dir=paste0(table_dir,"/result/",j,"/"),
      filter=j,
      padj.thresholds=padj_thresholds,
      abs.lfc.thresholds=abs_lfc_thresholds
    )
  }
}

############################################
#      Export lfcShrink Results
############################################

exportLfcShrinkResults = function(
  ddsMatrix,
  res,
  contrast,
  annotations,
  annotations.names,
  table.dir,
  filter=c("none","independent","standardized"),
  padj.thresholds=c(1),
  abs.lfc.thresholds=c(0)
) {

  #join by
  col = "ensembl_gene_id_version"

  #filter_suffix
  filter = match.arg(filter)
  filter_suffix = switch(filter,
    none="no-filt",
    independent="indep-filt",
    standardized="stndrd-filt"
  )

  #blind_suffix
  type_suffix = "apeglm"

  #result_lfcShrink
  coef=resultsNames(ddsMatrix)[length(resultsNames(ddsMatrix))]
  result_lfcShrink = suppressMessages(
    lfcShrink(
      dds= ddsMatrix,
      res = res,
      type = "apeglm",
      coef = coef,
      parallel = TRUE, 
      BPPARAM = SnowParam(workers=workers_available,progressbar=F),
      format = "DataFrame"                                    
    )
  )

  #lfcShrink
  result_lfcShrink = as.data.frame(result_lfcShrink) %>%
    rownames_to_column(col)

  #anno
  result_lfcShrink_anno=list()
  for (i in 1:length(annotations)) {
    if (!is.null(annotations[[i]])) {
      tmp = result_lfcShrink
      for (j in 1:length(annotations[[i]])) {
        tmp = tmp %>% left_join(annotations[[i]][[j]], by=col)
      }
      result_lfcShrink_anno[[i]] = tmp
    } else {
      result_lfcShrink_anno[[i]] = result_lfcShrink
    }
  }

  #filter result_lfcShrink for only significant genes
  result_lfcShrink_signif = list()
  filename = list()

  l = 1
  for(i in 1:length(result_lfcShrink_anno)) {
    for (j in 1:length(padj.thresholds)) {
      for (k in 1:length(abs.lfc.thresholds)) {
        result_lfcShrink_signif[[l]] = result_lfcShrink_anno[[i]] %>%
          dplyr::filter(padj < padj.thresholds[j]) %>%
          dplyr::filter(abs(log2FoldChange) > abs.lfc.thresholds[k])
        signif_suffix = paste0("padj",padj.thresholds[j],"_lfc",abs.lfc.thresholds[k])
        filename[[l]] = paste0("result-lfcShrink","_",filter_suffix,"_",annotations.names[i],"_",signif_suffix,".csv")
        l = l + 1
      }
    }
  }
  filename=unlist(filename)

  #make dir
  final.dir = paste0(table.dir,"/",contrast,"/")
  if(!dir.exists(final.dir)) { dir.create(final.dir,recursive = T) }

  #write tables
  for (i in 1:length(result_lfcShrink_signif)) {
    write.csv(x=result_lfcShrink_signif[[i]], file=paste0(final.dir,filename[i]), row.names = F)
  }
  print(contrast)
  print(filename)
}

#result_lfcShrink execute
for (i in 1:(length(contrasts)-1)) {
  for (j in filters) {
    exportLfcShrinkResults(
      ddsMatrix=ddsMatrices_complete[[j]][[i]],
      res=results_complete[[j]][[i]],
      contrast=contrasts[i],
      annotations=possible_annotations,
      annotations.names=names(possible_annotations),
      table.dir=paste0(table_dir,"/result_lfcShrink/",j,"/"),
      filter=j,
      padj.thresholds=padj_thresholds,
      abs.lfc.thresholds=abs_lfc_thresholds
    )
  }
}

############################################
#     Export Count Transformations
############################################

exportRegularizedLog = function(
  ddsMatrix,
  contrast,
  annotations,
  annotations.names,
  table.dir,
  filter=c("none","independent","standardized"),
  blind=TRUE
) {

  #join by
  col = "ensembl_gene_id_version"

  #filter_suffix
  filter = match.arg(filter)
  filter_suffix = switch(filter,
    none="no-filt",
    independent="indep-filt",
    standardized="stndrd-filt"
  )

  #blind_suffix
  blind_suffix = ifelse(blind, "blind", "not-blind")

  #rlog
  rlog = as.data.frame(assay(rlog(ddsMatrix, blind=blind))) %>%
    rownames_to_column(col)
  rlog_tidy = rlog %>%
    pivot_longer(!ensembl_gene_id_version,names_to="sample_id",values_to="counts")


  #anno
  rlog_anno=list()
  rlog_tidy_anno=list()

  for (i in 1:length(annotations)) {
    if (!is.null(annotations[[i]])) {
      rlog_tmp = rlog
      rlog_tidy_tmp = rlog_tidy
      for (j in 1:length(annotations[[i]])) {
        rlog_tmp = rlog_tmp %>% left_join(annotations[[i]][[j]], by=col)
        rlog_tidy_tmp = rlog_tidy_tmp %>% left_join(annotations[[i]][[j]], by=col)
      }
      rlog_anno[[i]] = rlog_tmp
      rlog_tidy_anno[[i]] = rlog_tidy_tmp
    } else {
      rlog_anno[[i]] = rlog
      rlog_tidy_anno[[i]] = rlog_tidy
    }
  }

  #filenames
  rlog_filename=list()
  rlog_tidy_filename=list()

  for(i in 1:length(rlog_anno)) {
    rlog_filename[[i]] = paste0("rlog","_",filter_suffix,"_",annotations.names[i],"_",blind_suffix,".csv")
    rlog_tidy_filename[[i]] = paste0("rlog-tidy","_",filter_suffix,"_",annotations.names[i],"_",blind_suffix,".csv")
  }

  rlog_filename = unlist(rlog_filename)
  rlog_tidy_filename = unlist(rlog_tidy_filename)

  #make dir
  final.dir = paste0(table.dir,"/",contrast,"/")
  if(!dir.exists(final.dir)) { dir.create(final.dir,recursive = T) }

  #write tables
  for (i in 1:length(rlog_anno)) {
    write.csv(x=rlog_anno[[i]], file=paste0(final.dir,rlog_filename[i]), row.names = F)
    write.csv(x=rlog_tidy_anno[[i]], file=paste0(final.dir,rlog_tidy_filename[i]), row.names = F)
  }
  print(contrast)
  print(rlog_filename)
  print(rlog_tidy_filename)
}

#rlog execute
for (i in 1:(length(contrasts)-1)) {
  for (j in blind_options) {
    for (k in filters) {
      exportRegularizedLog(
        ddsMatrix=ddsMatrices_complete[[k]][[i]],
        contrast=contrasts[i],
        annotations=possible_annotations,
        annotations.names=names(possible_annotations),
        table.dir=paste0(table_dir,"/rlog/",k,"/"),
        filter=k,
        blind=j
      )
    }
  }
}

#rlog execute : all
for (i in blind_options) {
  for (j in c("none", "standardized")) {
    exportRegularizedLog(
      ddsMatrix=ddsMatrices_all_complete[[j]],
      contrast="all",
      annotations=possible_annotations,
      annotations.names=names(possible_annotations),
      table.dir=paste0(table_dir,"/rlog/",j,"/"),
      filter=j,
      blind=i
    )
  }
}

exportVarianceStabilizingTransformation = function(
  ddsMatrix,
  contrast,
  annotations,
  annotations.names,
  table.dir,
  filter=c("none","independent","standardized"),
  blind=TRUE
) {

  #join by
  col = "ensembl_gene_id_version"

  #filter_suffix
  filter = match.arg(filter)
  filter_suffix = switch(filter,
    none="no-filt",
    independent="indep-filt",
    standardized="stndrd-filt"
  )

  #blind_suffix
  blind_suffix = ifelse(blind, "blind", "not-blind")

  #vst
  vst = as.data.frame(assay(vst(ddsMatrix, blind=blind))) %>%
    rownames_to_column(col)
  vst_tidy = vst %>%
    pivot_longer(!ensembl_gene_id_version,names_to="sample_id",values_to="counts")


  #anno
  vst_anno=list()
  vst_tidy_anno=list()

  for (i in 1:length(annotations)) {
    if (!is.null(annotations[[i]])) {
      vst_tmp = vst
      vst_tidy_tmp = vst_tidy
      for (j in 1:length(annotations[[i]])) {
        vst_tmp = vst_tmp %>% left_join(annotations[[i]][[j]], by=col)
        vst_tidy_tmp = vst_tidy_tmp %>% left_join(annotations[[i]][[j]], by=col)
      }
      vst_anno[[i]] = vst_tmp
      vst_tidy_anno[[i]] = vst_tidy_tmp
    } else {
      vst_anno[[i]] = vst
      vst_tidy_anno[[i]] = vst_tidy
    }
  }

  #filenames
  vst_filename=list()
  vst_tidy_filename=list()

  for(i in 1:length(vst_anno)) {
    vst_filename[[i]] = paste0("vst","_",filter_suffix,"_",annotations.names[i],"_",blind_suffix,".csv")
    vst_tidy_filename[[i]] = paste0("vst-tidy","_",filter_suffix,"_",annotations.names[i],"_",blind_suffix,".csv")
  }

  vst_filename = unlist(vst_filename)
  vst_tidy_filename = unlist(vst_tidy_filename)

  #make dir
  final.dir = paste0(table.dir,"/",contrast,"/")
  if(!dir.exists(final.dir)) { dir.create(final.dir,recursive = T) }

  #write tables
  for (i in 1:length(vst_anno)) {
    write.csv(x=vst_anno[[i]], file=paste0(final.dir,vst_filename[i]), row.names = F)
    write.csv(x=vst_tidy_anno[[i]], file=paste0(final.dir,vst_tidy_filename[i]), row.names = F)
  }
  print(contrast)
  print(vst_filename)
  print(vst_tidy_filename)
}

#vst execute
for (i in 1:(length(contrasts)-1)) {
  for (j in blind_options) {
    for (k in filters) {
      exportVarianceStabilizingTransformation(
        ddsMatrix=ddsMatrices_complete[[k]][[i]],
        contrast=contrasts[i],
        annotations=possible_annotations,
        annotations.names=names(possible_annotations),
        table.dir=paste0(table_dir,"/vst/",k,"/"),
        filter=k,
        blind=j
      )
    }
  }
}

#vst execute : all
for (i in blind_options) {
  for (j in c("none", "standardized")) {
    exportVarianceStabilizingTransformation(
      ddsMatrix=ddsMatrices_all_complete[[j]],
      contrast="all",
      annotations=possible_annotations,
      annotations.names=names(possible_annotations),
      table.dir=paste0(table_dir,"/vst/",j,"/"),
      filter=j,
      blind=i
    )
  }
}

exportFPKM = function(
  ddsMatrix,
  contrast,
  annotations,
  annotations.names,
  table.dir,
  filter=c("none","independent","standardized"),
  robust=TRUE,
  exons.by.gene
) {

  #join by
  col = "ensembl_gene_id_version"

  #filter_suffix
  filter = match.arg(filter)
  filter_suffix = switch(filter,
    none="no-filt",
    independent="indep-filt",
    standardized="stndrd-filt"
  )

  #robust_suffix
  robust_suffix = ifelse(robust, "robust", "not-robust")

  #fpkm
  ddsMatrix_FPKM = ddsMatrix
  genes_included = rownames(ddsMatrix)
  rowRanges(ddsMatrix_FPKM) = GRangesList(exons.by.gene[genes_included])
  fpkm = as.data.frame(fpkm(ddsMatrix_FPKM, robust=robust)) %>%
    rownames_to_column(col)
  fpkm_tidy = fpkm %>%
    pivot_longer(!ensembl_gene_id_version,names_to="sample_id",values_to="counts")

  #anno
  fpkm_anno=list()
  fpkm_tidy_anno=list()

  for (i in 1:length(annotations)) {
    if (!is.null(annotations[[i]])) {
      fpkm_tmp = fpkm
      fpkm_tidy_tmp = fpkm_tidy
      for (j in 1:length(annotations[[i]])) {
        fpkm_tmp = fpkm_tmp %>% left_join(annotations[[i]][[j]], by=col)
        fpkm_tidy_tmp = fpkm_tidy_tmp %>% left_join(annotations[[i]][[j]], by=col)
      }
      fpkm_anno[[i]] = fpkm_tmp
      fpkm_tidy_anno[[i]] = fpkm_tidy_tmp
    } else {
      fpkm_anno[[i]] = fpkm
      fpkm_tidy_anno[[i]] = fpkm_tidy
    }
  }

  #filenames
  fpkm_filename=list()
  fpkm_tidy_filename=list()

  for(i in 1:length(fpkm_anno)) {
    fpkm_filename[[i]] = paste0("fpkm","_",filter_suffix,"_",annotations.names[i],"_",robust_suffix,".csv")
    fpkm_tidy_filename[[i]] = paste0("fpkm-tidy","_",filter_suffix,"_",annotations.names[i],"_",robust_suffix,".csv")
  }

  fpkm_filename = unlist(fpkm_filename)
  fpkm_tidy_filename = unlist(fpkm_tidy_filename)

  #make dir
  final.dir = paste0(table.dir,"/",contrast,"/")
  if(!dir.exists(final.dir)) { dir.create(final.dir,recursive = T) }

  #write tables
  for (i in 1:length(fpkm_anno)) {
    write.csv(x=fpkm_anno[[i]], file=paste0(final.dir,fpkm_filename[i]), row.names = F)
    write.csv(x=fpkm_tidy_anno[[i]], file=paste0(final.dir,fpkm_tidy_filename[i]), row.names = F)
  }
  print(contrast)
  print(fpkm_filename)
  print(fpkm_tidy_filename)
}

#fpkm execute
for (i in 1:(length(contrasts)-1)) {
  for (j in robust_options) {
    for (k in filters) {
      exportFPKM(
        ddsMatrix=ddsMatrices_complete[[k]][[i]],
        contrast=contrasts[i],
        annotations=possible_annotations,
        annotations.names=names(possible_annotations),
        table.dir=paste0(table_dir,"/fpkm/",k,"/"),
        filter=k,
        robust=j,
        exons.by.gene = exons_by_gene
      )
    }
  }
}

#fpkm execute : all
for (i in robust_options) {
  for (j in c("none", "standardized")) {
    exportFPKM(
      ddsMatrix=ddsMatrices_all_complete[[j]],
      contrast="all",
      annotations=possible_annotations,
      annotations.names=names(possible_annotations),
      table.dir=paste0(table_dir,"/fpkm/",j,"/"),
      filter=j,
      robust=i,
      exons.by.gene = exons_by_gene
    )
  }
}

############################################
#      Export Counts
############################################


exportCounts = function(
  ddsMatrix,
  contrast,
  annotations,
  annotations.names,
  table.dir,
  filter=c("none","independent","standardized"),
  normalized=TRUE
) {

  #join by
  col = "ensembl_gene_id_version"

  #filter_suffix
  filter = match.arg(filter)
  filter_suffix = switch(filter,
    none="no-filt",
    independent="indep-filt",
    standardized="stndrd-filt"
  )

  #normalized_suffix
  normalized_suffix = ifelse(normalized, "size-norm", "not-norm")

  #count
  count = as.data.frame(counts(ddsMatrix, normalized=normalized)) %>%
    rownames_to_column(col)
  count_tidy = count %>%
    pivot_longer(!ensembl_gene_id_version,names_to="sample_id",values_to="counts")


  #anno
  count_anno=list()
  count_tidy_anno=list()

  for (i in 1:length(annotations)) {
    if (!is.null(annotations[[i]])) {
      count_tmp = count
      count_tidy_tmp = count_tidy
      for (j in 1:length(annotations[[i]])) {
        count_tmp = count_tmp %>% left_join(annotations[[i]][[j]], by=col)
        count_tidy_tmp = count_tidy_tmp %>% left_join(annotations[[i]][[j]], by=col)
      }
      count_anno[[i]] = count_tmp
      count_tidy_anno[[i]] = count_tidy_tmp
    } else {
      count_anno[[i]] = count
      count_tidy_anno[[i]] = count_tidy
    }
  }

  #filenames
  count_filename=list()
  count_tidy_filename=list()

  for(i in 1:length(count_anno)) {
    count_filename[[i]] = paste0("count","_",filter_suffix,"_",annotations.names[i],"_",normalized_suffix,".csv")
    count_tidy_filename[[i]] = paste0("count-tidy","_",filter_suffix,"_",annotations.names[i],"_",normalized_suffix,".csv")
  }

  count_filename = unlist(count_filename)
  count_tidy_filename = unlist(count_tidy_filename)

  #make dir
  final.dir = paste0(table.dir,"/",contrast,"/")
  if(!dir.exists(final.dir)) { dir.create(final.dir,recursive = T) }

  #write tables
  for (i in 1:length(count_anno)) {
    write.csv(x=count_anno[[i]], file=paste0(final.dir,count_filename[i]), row.names = F)
    write.csv(x=count_tidy_anno[[i]], file=paste0(final.dir,count_tidy_filename[i]), row.names = F)
  }
  print(contrast)
  print(count_filename)
  print(count_tidy_filename)
}

#count execute
for (i in 1:(length(contrasts)-1)) {
  for (j in normalized_options) {
    for (k in filters) {
      exportCounts(
        ddsMatrix=ddsMatrices_complete[[k]][[i]],
        contrast=contrasts[i],
        annotations=possible_annotations,
        annotations.names=names(possible_annotations),
        table.dir=paste0(table_dir,"/count/",k,"/"),
        filter=k,
        normalized=j
      )
    }
  }
}

#count execute : all
for (i in normalized_options) {
  for (j in c("none", "standardized")) {
    exportCounts(
      ddsMatrix=ddsMatrices_all_complete[[j]],
      contrast="all",
      annotations=possible_annotations,
      annotations.names=names(possible_annotations),
      table.dir=paste0(table_dir,"/count/",j,"/"),
      filter=j,
      normalized=i
    )
  }
}

############################################
#     Export Count Transformations
############################################

# exportZScore = function(
#   ddsMatrix,
#   contrast,
#   annotations,
#   annotations.names,
#   table.dir,
#   method=c("rlog","vst"),
#   filter=c("none","independent","standardized"),
#   blind=TRUE,
#   unique.groupings,
#   sampleLookup=sampleLookup
# ) {

#   #join by
#   col = "ensembl_gene_id_version"

#   #filter_suffix
#   filter=match.arg(filter)
#   filter_suffix = switch(filter,
#     none="no-filt",
#     independent="indep-filt",
#     standardized="stndrd-filt"
#   )

#   #method match
#   blind_suffix = ifelse(blind, "blind", "not-blind")

#   #zscore
#   method=match.arg(method)
#   zscore = switch(method,
#     rlog = t(scale(t(assay(rlog(ddsMatrix, blind=blind))))),
#     vst = t(scale(t(assay(vst(ddsMatrix, blind=blind)))))
#   )
#   zscore = as.data.frame(zscore) %>%
#     rownames_to_column(col)

#   #generate unique groups
#   y_list = list()
#   for (i in 1:length(unique.groupings)) {
#     selection=unique.groupings[[i]]
#     x = sampleLookup %>%
#       dplyr::select(selection) %>%
#       unique()
#     x$unique_groups = apply( x[, selection ] , 1 , paste , collapse = "_" )
    
#     y_list[[i]] = sampleLookup %>%
#       left_join(x)
#   }
#   sampleLookup_w_groups = bind_rows(y_list)

#   #unique_group averages
#   u_groups = levels(factor(sampleLookup_w_groups$unique_groups))
#   for(i in 1:length(u_groups)) {
#     group_ids = sampleLookup_w_groups$sample_id[sampleLookup_w_groups$unique_groups == u_groups[i]]

#     #with(data=sampleLookup_w_groups, expr=sample_id[])
#     zscore[[u_groups[i]]] = rowMeans(subset(zscore, select = group_ids), na.rm = TRUE)
#   }

#   zscore_tidy = zscore %>%
#     pivot_longer(!ensembl_gene_id_version,names_to="sample_id",values_to="zscore")


#   #anno
#   zscore_anno=list()
#   zscore_tidy_anno=list()

#   for (i in 1:length(annotations)) {
#     if (!is.null(annotations[[i]])) {
#       zscore_tmp = zscore
#       zscore_tidy_tmp = zscore_tidy
#       for (j in 1:length(annotations[[i]])) {
#         zscore_tmp = zscore_tmp %>% left_join(annotations[[i]][[j]], by=col)
#         zscore_tidy_tmp = zscore_tidy_tmp %>% left_join(annotations[[i]][[j]], by=col)
#       }
#       zscore_anno[[i]] = zscore_tmp
#       zscore_tidy_anno[[i]] = zscore_tidy_tmp
#     } else {
#       zscore_anno[[i]] = zscore
#       zscore_tidy_anno[[i]] = zscore_tidy
#     }
#   }

#   #filenames
#   zscore_filename=list()
#   zscore_tidy_filename=list()

#   for(i in 1:length(zscore_anno)) {
#     zscore_filename[[i]] = paste0("zscore","_",method,"_",filter_suffix,"_",annotations.names[i],"_",blind_suffix,".csv")
#     zscore_tidy_filename[[i]] = paste0("zscore-tidy","_",method,"_",filter_suffix,"_",annotations.names[i],"_",blind_suffix,".csv")
#   }

#   zscore_filename = unlist(zscore_filename)
#   zscore_tidy_filename = unlist(zscore_tidy_filename)

#   #make dir
#   final.dir = paste0(table.dir,"/",contrast,"/")
#   if(!dir.exists(final.dir)) { dir.create(final.dir,recursive = T) }

#   #write tables
#   for (i in 1:length(zscore_anno)) {
#     write.csv(x=zscore_anno[[i]], file=paste0(final.dir,zscore_filename[i]), row.names = F)
#     write.csv(x=zscore_tidy_anno[[i]], file=paste0(final.dir,zscore_tidy_filename[i]), row.names = F)
#   }
#   print(contrast)
#   print(zscore_filename)
#   print(zscore_tidy_filename)
# }

# #FIX_HARDCODE
# unique_groupings = list(c("tissue","treatment"),c("tissue","treatment_simple"),c("treatment"),c("treatment_simple"))

# #count execute
# for (i in 1:(length(contrasts)-1)) {
#   for (j in blind_options) {
#     for (k in filters) {
#       for (l in c("vst","rlog")) {
#         exportZScore(
#           ddsMatrix=ddsMatrices_complete[[k]][[i]],
#           contrast=contrasts[i],
#           annotations=possible_annotations,
#           annotations.names=names(possible_annotations),
#           table.dir=paste0(table_dir,"/zscore/",k,"/"),
#           filter=k,
#           blind=j,
#           method=l,
#           unique.groupings = unique_groupings
#         )
#       }
#     }
#   }
# }

# #count execute : all
# for (i in c(FALSE)) {
#   for (j in c("none", "standardized")) {
#     for (l in c("rlog")) {
#       exportZScore(
#         ddsMatrix=ddsMatrices_all_complete[[j]],
#         contrast="all",
#         annotations=possible_annotations,
#         annotations.names=names(possible_annotations),
#         table.dir=paste0(table_dir,"/zscore/",j,"/"),
#         filter=j,
#         blind=i,
#         unique.groupings = unique_groupings,
#         method=l,
#         sampleLookup=sampleLookup
#       )
#     }
#   }
# }

