#!/usr/bin/env Rscript

###---PACKAGES---###

if (!require("pacman")) { install.packages("pacman", repos='http://cran.us.r-project.org') }
library(pacman)

#required packages
required_packages = c(
  "tidyverse",
  "grid",
  "ggplotify",
  "svglite"
)

github_packages = c(
  "slowkow/ggrepel"
)

#load packages
pacman::p_load(
  char=required_packages,
  install=TRUE,
  character.only=TRUE,
  try.bioconductor=TRUE,
  update.bioconductor=TRUE
)

#load github packages
pacman::p_load_gh(
  char = github_packages,
  update = getOption("pac_update"),
  dependencies = TRUE
)

###---GLOBALS---###

genesToLabel=NULL
#if you want to label points, just make a list
#i.e. c("MITF")

padj_threshold = 0.01
abs_lfc_threshold = 1

###---CLEAN DATA---###

data.in = "../data/in/dge/featureCounts_deseq2/table/result_lfcShrink/standardized/"
data.out = "../data/out/"
plot.out = "../plot/"

#contrasts
contrasts = basename(list.dirs(path = data.in))[2:length(list.dirs(path = data.in))]

#results
results=list()
for (i in 1:length(contrasts)) {
	results[[i]] = read.csv(file = paste0(data.in,contrasts[i],"/result-lfcShrink_stndrd-filt_anno-basic_padj1_lfc0.csv"), header = TRUE)
}
names(results) = contrasts

clean_data = function(result,genesToLabel = NULL) {

	result.no.zeroes = result %>%
		dplyr::filter(padj != 0)

	min.padj.value = min(result.no.zeroes$padj)

    ret=result %>%
      dplyr::mutate(padj=ifelse(padj == 0, min.padj.value, padj)) %>%
      dplyr::mutate(color_group="Not Significant") %>%
      dplyr::mutate(color_group=ifelse(padj < padj_threshold & log2FoldChange < -abs_lfc_threshold, "Down-Regulated Genes", color_group)) %>%
      dplyr::mutate(color_group=ifelse(padj < padj_threshold & log2FoldChange > abs_lfc_threshold, "Up-Regulated Genes", color_group))

    if (!is.null(genesToLabel)) {
    ret = ret %>%	
      dplyr::mutate(label=ifelse(external_gene_name %in% genesToLabel, external_gene_name, NA))
    }
  return(ret)
}

cleaned_inputs = list()
for (i in 1:length(results)) {
	cleaned_inputs[[i]] = clean_data(result=results[[i]],genesToLabel = genesToLabel)
}
names(cleaned_inputs) = contrasts

###---PLOT---###

x_lab=expression(-log[2](Fold~Change))
y_lab=expression(-log[10](italic(q)-value))

volcano_plot = function(cleaned_input,x_lab=x_lab,y_lab=y_lab) {
  cols=c("Down-Regulated Genes" = "#234463","Up-Regulated Genes" = "#781e1e", "Not Significant" = "gray50")
  cols2=c("Down-Regulated Genes" = "#f1f8ff", "Up-Regulated Genes" = "#fff6f6")

  max_lfc=max(abs(cleaned_input$log2FoldChange))
  max_padj=max(-log10(cleaned_input$padj),na.rm = TRUE)
  
  plot=ggplot(data=cleaned_input, mapping=aes(x=log2FoldChange,y=-log10(padj))) +
    #down
    geom_rect(
      fill = cols2[1],
      xmin = -Inf, 
      xmax = -abs_lfc_threshold,
      ymin = -log10(padj_threshold),
      ymax = Inf
    ) +
    #up
    geom_rect(
        fill = cols2[2],
        xmin = abs_lfc_threshold, 
        xmax = Inf,
        ymin = -log10(padj_threshold),
        ymax = Inf
      ) +
    geom_point(mapping=aes(color=color_group,fill=color_group),alpha=0.5) +
    geom_hline(yintercept=-log10(padj_threshold), color='black', size=0.5, linetype = "dashed") +
    geom_vline(xintercept=abs_lfc_threshold, color='black', size=0.5, linetype = "dashed") +
    geom_vline(xintercept=-abs_lfc_threshold, color='black', size=0.5, linetype = "dashed") +
    scale_color_manual(values=cols,guide=FALSE) +
    scale_fill_manual(values=cols,guide=FALSE) +
    scale_x_continuous(limits = c(-max_lfc-1,max_lfc+1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0,(max_padj+1)), expand = c(0, 0)) +
    labs(x=x_lab,y=y_lab) + {
    	if ("label" %in% names(cleaned_input)) {
	    	ggrepel::geom_label_repel(
		      min.segment.length = 0,
		      mapping = aes(label=label)
	    	) 
    	}
    } +
    theme_classic() +
    theme(
      axis.title=element_text(size=12),
      strip.text=element_text(size=12, color = "white", face="bold"),
      axis.text=element_text(size=12),
      axis.line = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1.)
    )
  #aes(fill = as.factor(geneset))
  
  #strip colors
  # g = ggplot_gtable(ggplot_build(plot))
  # striprt = which( grepl('strip-t', g$layout$name) )
  # fills = c("#234463","#781e1e")
  # k = 1
  # for (i in striprt) {
  #   j = which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  #   g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill = fills[k]
  #   k = k+1
  # }
  # return(ggplotify::as.ggplot(g))
  return(plot)
}

plots = list()
for (i in 1:length(results)) {
	plots[[i]] = volcano_plot(cleaned_input = cleaned_inputs[[i]],x_lab = x_lab,y_lab = y_lab)
}
names(plots) = contrasts

for(i in 1:length(plots)) {
	ggsave(filename=paste0(plot.out,"volcano_",names(plots)[i],".png"),plot=plots[[i]],device="png",dpi=320,width=6,height=6)
	ggsave(filename=paste0(plot.out,"volcano_",names(plots)[i],".svg"),plot=plots[[i]],device="svg",dpi=320,width=6,height=6)
	ggsave(filename=paste0(plot.out,"volcano_",names(plots)[i],".pdf"),plot=plots[[i]],device="pdf",dpi=320,width=6,height=6)
}