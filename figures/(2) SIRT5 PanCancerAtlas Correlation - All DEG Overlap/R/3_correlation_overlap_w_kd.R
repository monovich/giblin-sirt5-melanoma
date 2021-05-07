#!/usr/bin/env Rscript

###---PACKAGES---###

if (!require("pacman")) { install.packages("pacman", repos='http://cran.us.r-project.org') }
library(pacman)

#required packages
required_packages = c(
  "tidyverse",
  "grid",
  "ggplotify",
  "svglite",
  "AnnotationDbi",
  "org.Hs.eg.db"
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
padj_rna_threshold = 0.01
lfc_rna_threshold = 1
padj_rna_label_threshold = 1
lfc_rna_label_threshold = 0
padj_rho_threshold = 0.01

data.in = "../data/in/"
data.in.long = "../data/in/dge/featureCounts_deseq2/table/result_lfcShrink/standardized/sirt5_kd_over_sirt5_nt/"
data.out = "../data/out/"
plot.out = "../plot/"

###---IMPORT DATA---#

#correlation data
correlation_w_SIRT5 = read.csv(file = paste0(data.out,"correlation_w_SIRT5.csv")) %>%
	as_tibble() %>%
	dplyr::mutate(fontface = "italic")
correlation_w_SIRT5 = correlation_w_SIRT5[
!is.na(correlation_w_SIRT5$spearman.rho) & 
!is.na(correlation_w_SIRT5$spearman.FDR) &
!is.na(correlation_w_SIRT5$pearson.cor) &
!is.na(correlation_w_SIRT5$pearson.FDR)
,]

correlation_w_SIRT5$ENTREZID = as.character(correlation_w_SIRT5$ENTREZID)

#non-sirt5 results filtered
result = read.csv(file = paste0(data.in.long,"result-lfcShrink_stndrd-filt_anno-basic_padj1_lfc0.csv")) %>%
	as_tibble() %>%
  dplyr::filter(external_gene_name != "SIRT5") %>%
	dplyr::filter(padj < padj_rna_threshold & abs(log2FoldChange) > lfc_rna_threshold)

###---ENTREZ MAPPINGS---###

# Return the Entrez mapping for a set of genes
annotations_orgDb = AnnotationDbi::select(
    org.Hs.eg.db,
    keys = result$ensembl_gene_id,
    columns = c("SYMBOL", "ENTREZID","GENENAME"),
    keytype = "ENSEMBL"
  ) %>%
  as_tibble()

annotations_orgDb_original = annotations_orgDb

# Determine the indices for the non-NA genes
non_na_idx = which(is.na(annotations_orgDb$SYMBOL) == FALSE)

# Return only the genes with annotations using indices
annotations_orgDb = annotations_orgDb[non_na_idx, ]

# Determine the indices for the non-duplicated genes
non_duplicates_idx <- which(duplicated(annotations_orgDb$SYMBOL) == FALSE)

# Return only the non-duplicated genes using indices
annotations_orgDb <- annotations_orgDb[non_duplicates_idx, ]

# Determine the indices for the non-duplicated genes
non_duplicates_idx <- which(duplicated(annotations_orgDb$ENTREZID) == FALSE)

# Return only the non-duplicated genes using indices
annotations_orgDb <- annotations_orgDb[non_duplicates_idx, ]

# Determine the indices for the non-duplicated genes
non_duplicates_idx <- which(duplicated(annotations_orgDb$GENENAME) == FALSE)

# Return only the non-duplicated genes using indices
annotations_orgDb <- annotations_orgDb[non_duplicates_idx, ]

# Determine the indices for the non-duplicated genes
non_duplicates_idx <- which(duplicated(annotations_orgDb$ENSEMBL) == FALSE)

# Return only the non-duplicated genes using indices
annotations_orgDb <- annotations_orgDb[non_duplicates_idx, ]

annotations_orgDb = annotations_orgDb %>% setNames(c("ensembl_gene_id","SYMBOL","ENTREZID","GENENAME"))

result_w_entrez = result %>% 
  left_join(annotations_orgDb) %>%
  dplyr::filter(!is.na(ENTREZID))

#oncogenes
oncogenes_data = read.table(file = paste0(data.in,"ongene_human.txt"), sep = '\t', header = TRUE, quote = "")
oncogenes_entrezgene_ids = oncogenes_data$OncogeneID

#melanoma subtype signatures
melanoma_subtype_signatures = read.csv(file = paste0(data.in,"subtype_signatures_updated.csv"))
melanoma_external_gene_names = melanoma_subtype_signatures$Gene

###---PRE-PROCESS DATA---###

#subset correlation by significant SIRT5 correlations that match the direction of change in KD
# i.e. if down in KD, should have a positive correlation (SIRT5 is also down)
# i.e. if up in KD, should have a negative correlation (SIRT5 is down in KD, so needs to be) 
correlation_by_result = result_w_entrez %>%
	left_join(correlation_w_SIRT5) %>%
	dplyr::filter((spearman.rho > 0 & log2FoldChange < 0) | (spearman.rho < 0 & log2FoldChange > 0)) %>%
	dplyr::mutate(fontface = ifelse(external_gene_name %in% melanoma_external_gene_names, "bold.italic",fontface))


#filter to get genes we want to label
cbr_filtered = correlation_by_result %>%
	dplyr::filter(spearman.FDR < padj_rho_threshold) %>%
	dplyr::filter(
		ENTREZID %in% oncogenes_entrezgene_ids | 
		external_gene_name %in% melanoma_external_gene_names |
		abs(spearman.rho) > 0.35
	)

cbr_filtered_labels_to_omit = cbr_filtered %>%
  dplyr::filter(abs(spearman.rho) <= 0.35 & padj >= padj_rna_label_threshold & abs(log2FoldChange) <= lfc_rna_label_threshold)
ids_omit = cbr_filtered_labels_to_omit$ENTREZID

cbr_filtered = cbr_filtered %>%
  dplyr::filter(!ENTREZID %in% ids_omit)

entrezgene_ids_to_label = cbr_filtered$ENTREZID


###---FUNCTIONS---###

clean_data = function(result,entrezgene_label=c(),entrezgene_filter=NULL) {
	if (is.null(entrezgene_filter)) {
		ret = result
	} else {
		ret = result %>%
			dplyr::filter(ENTREZID %in% entrezgene_filter)
	}
    ret=ret %>%
      dplyr::mutate(color_group="Not Significant") %>%
      dplyr::mutate(color_group=ifelse(spearman.FDR < padj_rho_threshold & spearman.rho < 0 , "Negatively Correlated Genes", color_group)) %>%
      dplyr::mutate(color_group=ifelse(spearman.FDR < padj_rho_threshold & spearman.rho > 0, "Positively Correlated Genes", color_group)) %>%
      dplyr::mutate(left_even_label=ifelse(
      	(ENTREZID %in% entrezgene_label) & 
      	(spearman.rho < 0) & 
      	(spearman.FDR < padj_rho_threshold) &
      	(row_number() %% 2 == 0), external_gene_name, NA)) %>%
      dplyr::mutate(left_odd_label=ifelse(
      	(ENTREZID %in% entrezgene_label) & 
      	(spearman.rho < 0) & 
      	(spearman.FDR < padj_rho_threshold) &
      	(row_number() %% 2 == 1), external_gene_name, NA)) %>%
      dplyr::mutate(right_even_label=ifelse(
      	(ENTREZID %in% entrezgene_label) & 
      	(spearman.rho > 0) & 
      	(spearman.FDR < padj_rho_threshold) &
      	(row_number() %% 2 == 0), external_gene_name, NA)) %>%
      dplyr::mutate(right_odd_label=ifelse(
      	(ENTREZID %in% entrezgene_label) & 
      	(spearman.rho > 0) & 
      	(spearman.FDR < padj_rho_threshold) &
      	(row_number() %% 2 == 1), external_gene_name, NA))
  return(ret)
}

full_plot_data = clean_data(result=correlation_w_SIRT5)
subset_plot_data = clean_data(result=correlation_by_result,entrezgene_label=entrezgene_ids_to_label)

###---PLOT---###

x_lab=expression("Spearman's"~italic(rho))
y_lab=expression(-log[10](italic(q)-value))

set.seed(42)
correlation_volcano_plot = function(cleaned_input,x_lab,y_lab) {
  cols=c("Negatively Correlated Genes" = "#234463","Positively Correlated Genes" = "#781e1e", "Not Significant" = "gray50")
  cols2=c("Negatively Correlated Genes" = "#f1f8ff", "Positively Correlated Genes" = "#fff6f6")
  max_padj=max(-log10(cleaned_input$spearman.FDR),na.rm = TRUE)
  
  plot=ggplot(data=cleaned_input %>% arrange(-spearman.FDR), mapping=aes(x=spearman.rho,y=-log10(spearman.FDR))) +
    #down
    geom_rect(
      fill = cols2[1],
      xmin = -Inf, 
      xmax = 0,
      ymin = -log10(padj_rho_threshold),
      ymax = Inf
    ) +
    #up
    geom_rect(
        fill = cols2[2],
        xmin = 0, 
        xmax = Inf,
        ymin = -log10(padj_rho_threshold),
        ymax = Inf
    ) +
    geom_point(mapping=aes(color=color_group,fill=color_group),alpha=0.5) +
    geom_hline(yintercept=-log10(padj_rho_threshold), color='black', size=0.5, linetype = "dashed") +
    geom_vline(xintercept=0, color='black', size=0.5, linetype = "dashed") +
    scale_color_manual(values=cols,guide=FALSE) +
    scale_fill_manual(values=cols,guide=FALSE) +
    scale_x_continuous(limits = c(-0.75,0.75), breaks = c(-0.6,-0.4,-0.2,0,0.2,0.4,0.6), expand = c(0,0)) +
    scale_y_continuous(limits = c(0,(max_padj+1)), expand = c(0, 0)) +
    labs(x=x_lab,y=y_lab) +
    #right even (left) labels
    ggrepel::geom_text_repel(
      data=cleaned_input %>% arrange(-spearman.FDR),
      xlim  = c(0,NA),
      ylim = c(-log10(padj_rho_threshold-0.005),NA),
      nudge_x = -0.15,
      hjust=1,
      min.segment.length = 0,
      segment.square  = TRUE,
      segment.inflect = TRUE,
      segment.curvature = -1e-20,
      segment.ncp = 3,
      direction = "y",
      mapping = aes(label = right_even_label, fontface = fontface),
      box.padding = unit(0.1, "lines"),
      point.padding = unit(0.3, "lines"),
      size = 2,
      max.iter = 1e7,
      max.time = 2
    ) +
    #right odd (right) labels
    ggrepel::geom_text_repel(
      #xlim  = c(0.1,NA),
      data=cleaned_input %>% arrange(-spearman.FDR),
      ylim = c(-log10(padj_rho_threshold-0.005),NA),
      nudge_x = 0.15,
      hjust=0,
      min.segment.length = 0,
      segment.square  = TRUE,
      segment.inflect = TRUE,
      segment.curvature = -1e-20,
      segment.ncp = 3,
      direction = "y",
      mapping = aes(label = right_odd_label, fontface = fontface),
      box.padding = unit(0.1, "lines"),
      point.padding = unit(0.3, "lines"),
      size = 2,
      max.iter = 1e7,
      max.time = 2
    ) +
    #left even (left) labels
    ggrepel::geom_text_repel(
      #xlim  = c(NA,-0.05),
      data=cleaned_input %>% arrange(-spearman.FDR),
      ylim = c(-log10(padj_rho_threshold-0.005),NA),
      nudge_x = -0.15,
      hjust=1,
      min.segment.length = 0,
      segment.square  = TRUE,
      segment.inflect = TRUE,
      segment.curvature = -1e-20,
      segment.ncp = 3,
      direction = "y",
      mapping = aes(label = left_even_label, fontface = fontface),
      box.padding = unit(0.1, "lines"),
      point.padding = unit(0.3, "lines"),
      size = 2,      
      max.iter = 1e7,
      max.time = 2
    ) +
    #left odd (right) labels
    ggrepel::geom_text_repel(
      data=cleaned_input %>% arrange(-spearman.FDR),
      xlim  = c(NA,0),
      ylim = c(-log10(padj_rho_threshold-0.005),NA),
      nudge_x = 0.15,
      hjust=0,
      min.segment.length = 0,
      segment.square  = TRUE,
      segment.inflect = TRUE,
      segment.curvature = -1e-20,
      segment.ncp = 3,
      direction = "y",
      mapping = aes(label = left_odd_label, fontface = fontface),
      box.padding = unit(0.1, "lines"),
      point.padding = unit(0.3, "lines"),
      size = 2,      
      max.iter = 1e7,
      max.time = 2
    ) +
    annotation_custom(
      grob = textGrob(label = expression(bold(Negatively~Correlated~DEGs)), hjust = 0.5, gp = gpar(cex = 1)),
      ymin = max_padj + 1.6,
      ymax = max_padj + 1.6,
      xmin = -0.375,    
      xmax = -0.375
    ) +
    annotation_custom(
      grob = textGrob(label = expression(bold(Positively~Correlated~DEGs)), hjust = 0.5, gp = gpar(cex = 1)),
      ymin = max_padj + 1.6,
      ymax = max_padj + 1.6,
      xmin = 0.375,    
      xmax = 0.375
    ) +
    theme_classic(
    ) +
    theme(
      axis.title=element_text(size=12),
      strip.text=element_text(size=12, color = "white", face="bold"),
      axis.text=element_text(size=12),
      axis.line = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      aspect.ratio=1
    )

    #turn off boundary clipping
  gt <- ggplot_gtable(ggplot_build(plot))
	gt$layout$clip[gt$layout$name == "panel"] <- "off"
	return(gt)
}

full_plot = correlation_volcano_plot(cleaned_input=full_plot_data,x_lab=x_lab,y_lab=y_lab)
name="SIRT5-full-pancancer-correlation-volcano-plot"
ggsave(filename=paste0(plot.out,name,".png"),plot=full_plot,device="png",dpi=320,width=6.5,height=7)
ggsave(filename=paste0(plot.out,name,".svg"),plot=full_plot,device="svg",dpi=320,width=6.5,height=7)
ggsave(filename=paste0(plot.out,name,".pdf"),plot=full_plot,device="pdf",dpi=320,width=6.5,height=7)


subset_plot = correlation_volcano_plot(cleaned_input=subset_plot_data,x_lab=x_lab,y_lab=y_lab)
name="SIRT5-subset-pancancer-correlation-volcano-plot"
ggsave(filename=paste0(plot.out,name,".png"),plot=subset_plot,device="png",dpi=320,width=6.5,height=7)
ggsave(filename=paste0(plot.out,name,".svg"),plot=subset_plot,device="svg",dpi=320,width=6.5,height=7)
ggsave(filename=paste0(plot.out,name,".pdf"),plot=subset_plot,device="pdf",dpi=320,width=6.5,height=7)