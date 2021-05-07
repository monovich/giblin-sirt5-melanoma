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

###---GLOBAL CONFIG---###
ih_pvalue_threshold = 0.01
padj_threshold = 0.01
lfc_rna_threshold = 1

###---FUNCTIONS---###

data.in = "../data/in/"
data.in.long = "../data/in/dge/featureCounts_deseq2/table/result_lfcShrink/standardized/sirt5_kd_over_sirt5_nt/"
data.out = "../data/out/"
plot = "../plot/"

down_input=paste0(data.out,"down_genes_bart_results.txt")
up_input=paste0(data.out,"up_genes_bart_results.txt")

result_filtered = read.csv(paste0(data.in.long,"result-lfcShrink_stndrd-filt_anno-basic_padj1_lfc0.csv")) %>%
  dplyr::filter(external_gene_name != "SIRT5") %>%
  dplyr::filter(padj < padj_threshold) %>%
  dplyr::filter(abs(log2FoldChange) > lfc_rna_threshold)

#oncogenes
oncogenes_data = read.table(file = paste0(data.in,"ongene_human.txt"), sep = '\t', header = TRUE) %>%
    as_tibble()
oncogenes_names = oncogenes_data$OncogeneName

#melanoma subtype signatures
melanoma_subtype_signatures = read.csv(file = paste0(data.in,"subtype_signatures_updated.csv")) %>%
    as_tibble()
melanoma_external_gene_names = melanoma_subtype_signatures$Gene

genes=result_filtered$external_gene_name
genes=c(genes,melanoma_external_gene_names,oncogenes_names)

genes_bold=c(result_filtered$external_gene_name)

set.seed(42)

clean_data = function(up_input=NULL,down_input=NULL,genes=c(),genes.bold=c()) {
    up=NULL
    down=NULL
    ret=NULL
    if(!is.null(up_input)) {
        up=read.table(header=TRUE,file=up_input,sep='\t') %>%
            dplyr::mutate(geneset="Factors Predicted From Up-Regulated Genes") %>%
            dplyr::mutate(color_group=ifelse(irwin_hall_pvalue < ih_pvalue_threshold, "Factors Predicted From Up-Regulated Genes", "Not Significant")) %>%
            dplyr::mutate(label=ifelse((TF %in% genes) & (irwin_hall_pvalue < ih_pvalue_threshold), TF, NA)) %>%
            dplyr::mutate(fontface = ifelse((TF %in% genes_bold), "bold.italic","italic")) %>%
            dplyr::filter(zscore > 0)
    }
    if(!is.null(down_input)) {
        down=read.table(header=TRUE,file=down_input,sep='\t') %>%
            dplyr::mutate(geneset="Factors Predicted From Down-Regulated Genes") %>%
            dplyr::mutate(color_group=ifelse(irwin_hall_pvalue < ih_pvalue_threshold, "Factors Predicted From Down-Regulated Genes", "Not Significant")) %>%
            dplyr::mutate(label=ifelse((TF %in% genes) & (irwin_hall_pvalue < ih_pvalue_threshold), TF, NA)) %>%
            dplyr::mutate(fontface = ifelse((TF %in% genes_bold), "bold.italic","italic")) %>%
            dplyr::filter(zscore > 0)
    }
    if(!is.null(up) & is.null(down)) {
        #UP ONLY
        ret=up %>% as_tibble()
    } else if (is.null(up) & !is.null(down)) {
        #DOWN ONLY
        ret=down %>% as_tibble()
    } else if (!is.null(up) & !is.null(down)) {
        ret=rbind(up,down) %>% as_tibble()
    }
    return(ret)
}

x = clean_data(up_input = up_input, down_input = down_input,genes=genes,genes.bold=genes_bold)

signif_ovr_effect = function(cleaned_input) {
    expression1=expression(italic(z)-score)
    expression2=expression(-log[10](Irwin~Hall~italic(p)-value))
    cols=c("Factors Predicted From Down-Regulated Genes" = "#234463","Factors Predicted From Up-Regulated Genes" = "#781e1e", "Not Significant" = "gray50")
    cols2=c("Factors Predicted From Down-Regulated Genes" = "#f1f8ff", "Factors Predicted From Up-Regulated Genes" = "#fff6f6", "Not Significant" = "gray50")
    df = cleaned_input %>% dplyr::filter(-log10(irwin_hall_pvalue) > 1.85)
    max_zscore=max(df$zscore)
    median_zscore=median(df$zscore)
    max_ih=max(-log10(df$irwin_hall_pvalue))
    tmpplot=ggplot(data=cleaned_input, mapping=aes(x=zscore,y=-log10(irwin_hall_pvalue))) +
        geom_rect(
            mapping=aes(fill = geneset),
            xmin = -Inf, 
            xmax = Inf,
            ymin = 2,
            ymax = Inf
        ) +
        geom_point(mapping=aes(color=color_group,fill=color_group),alpha=0.5) +
        geom_hline(yintercept=range(-log10(0.01)), color='black', size=0.5, linetype = "dashed") +
        scale_color_manual(values=cols,guide=FALSE) +
        scale_fill_manual(values=cols2,guide=FALSE) +
        scale_x_continuous(limits = c(0,max_zscore+0.25), expand = c(0, 0)) +
        scale_y_continuous(limits = c(0,(max_ih+1)), expand = c(0, 0)) +
        labs(x=expression1,y=expression2) +
        ggrepel::geom_label_repel(
          nudge_x = -0.3,
          ylim = c(-log10(0.01),NA),
          hjust=0.5,
          min.segment.length = 0,
          segment.square  = TRUE,
          segment.inflect = TRUE,
          segment.curvature = -1e-20,
          segment.ncp = 3,
          fill=alpha("white",0.85),
          mapping = aes(label = label, fontface = fontface),
          box.padding = unit(0.1, "lines"),
          point.padding = unit(0.3, "lines"),
          size = 2,      
          max.iter = 1e7,
          max.time = 2
        ) +
        theme_classic() +
        theme(
            axis.title=element_text(size=12),
            strip.text=element_text(size=12, color = "white", face="bold"),
            axis.text=element_text(size=12),
            axis.line = element_blank(),
            panel.border = element_rect(color = "black", fill = NA, size = 1.)
        ) +
        facet_grid(cols = vars(geneset))

    #strip colors
    g = ggplot_gtable(ggplot_build(tmpplot))
    striprt = which( grepl('strip-t', g$layout$name) )
    fills = c("#234463","#781e1e")
    k = 1
    for (i in striprt) {
        j = which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
        g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill = fills[k]
        k = k+1
    }
    return(g)
}

bart_plot = signif_ovr_effect(cleaned_input=x)

#plot

name="BART_sirt5_kd_over_sirt5_nt"
ggsave(filename=paste0(plot,name,".png"),plot=bart_plot,device="png",dpi=320,width=10,height=7)
ggsave(filename=paste0(plot,name,".svg"),plot=bart_plot,device="svg",dpi=320,width=10,height=7)
ggsave(filename=paste0(plot,name,".pdf"),plot=bart_plot,device="pdf",dpi=320,width=10,height=7)
