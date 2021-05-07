#!/usr/bin/env Rscript

###---PACKAGES---###

if (!require("pacman")) { install.packages("pacman", repos='http://cran.us.r-project.org') }
library(pacman)

#required packages
required_packages = c(
  "tidyverse",
  "egg"
)

#load packages
pacman::p_load(
  char=required_packages,
  install=TRUE,
  character.only=TRUE,
  try.bioconductor=TRUE,
  update.bioconductor=TRUE
)

###---GLOBALS---###

set.seed(42)
padj_threshold = 0.01
lfc_rna_threshold = 1

###---FUNCTIONS---###

data.in = "../data/in/"
data.in.result = "../data/in/dge/featureCounts_deseq2/table/result_lfcShrink/standardized/sirt5_kd_over_sirt5_nt/"
data.in.fpkm = "../data/in/dge/featureCounts_deseq2/table/fpkm/standardized/all/"
data.out = "../data/out/"
plot = "../plot/"

#result
result = read.csv(file = paste0(data.in.result,"result-lfcShrink_stndrd-filt_anno-basic_padj1_lfc0.csv")) %>%
	dplyr::filter(padj < padj_threshold) %>%
  dplyr::filter(abs(log2FoldChange) > lfc_rna_threshold)

#sampleLookup
sampleLookup = read.csv(file = paste0(data.in,"dge/sampleLookup.csv")) 

#genes to plot
genes = result$external_gene_name

#fpkm
fpkm = read.csv(file = paste0(data.in.fpkm,"fpkm-tidy_stndrd-filt_anno-basic_robust.csv"))

for(k in 1:length(genes)) {
  gene = genes[k]
	fpkm_filterd = fpkm %>%
		dplyr::filter(external_gene_name == gene) %>%
		left_join(sampleLookup)
	plots = list()
	tissues = levels(factor(fpkm_filterd$tissue))
	treatment_labels = c("C","KD1","KD2")
	treatments = c("sirt5_nt","sirt5_kd_546","sirt5_kd_547")
	alpha=0.5
	colors = c("sirt5_nt"=alpha("blue",alpha),"sirt5_kd_546"=alpha("red",alpha),"sirt5_kd_547"=alpha("yellow",alpha))

	for (i in 1:length(tissues)) {
		dat = fpkm_filterd %>% dplyr::filter(tissue == tissues[i])
		annotations=data.frame(
						xpos=c(0.7),
						ypos=c(max(dat$counts)*1.3),
						annotateText = c(gene)
					)
		plots[[i]] = ggplot(data=dat,
							mapping = aes(x=factor(treatment,levels=treatments),y=counts)) +
			stat_boxplot(geom = "errorbar", width = 0.5) +
			geom_boxplot(outlier.shape = NA) +
			geom_jitter(aes(fill=treatment),color="black",size=3,shape=21,stroke=1) +
			scale_x_discrete(breaks=treatments, labels=treatment_labels) +
			#scale_color_manual(values=colors) +
			scale_fill_manual(values=colors,breaks=treatments,labels=treatment_labels) +
			#y=expression(log[10](FPKM))
			labs(x="",y="FPKM",title=tissues[i]) +
			lims(y=c(NA,max(dat$counts)*1.3)) +
			theme_classic() +
			theme(
	            axis.title=element_text(size=18,face="bold"),
	            axis.text=element_text(size=14),
	            axis.line = element_blank(),
	            panel.border = element_rect(color = "black", fill = NA, size = 1.),
	            title=element_text(size=20,face="italic"),
	            legend.text=element_text(size=14),
	            legend.title=element_blank()
			) +
			geom_text(
				data=annotations,
				mapping=aes(
					x=xpos,
					y=ypos,
					label=annotateText
				),
				hjust="inward",
				vjust="inward",
				size = 6,
				fontface="bold"
			)

		extension="svg"
		ggsave(
			filename=paste0(plot,tissues[i],"_",gene,".",extension),
			device=extension,
			plot=plots[[i]],
			width=4,
			height=7
		)

		extension="png"
		ggsave(
			filename=paste0(plot,tissues[i],"_",gene,".",extension),
			device=extension,
			plot=plots[[i]],
			width=4,
			height=7
		)
	}

	plots_fixed = list()
	plots_fixed[[1]] = plots[[1]] + theme(legend.position = "none")
	plots_fixed[[2]] = plots[[2]] + theme(axis.title.y=element_blank(),legend.position = "none")
	plots_fixed[[3]] = plots[[3]] + theme(axis.title.y=element_blank())

	extension="svg"
	arrangment = ggarrange(plots=plots_fixed,widths=c(4,4,4),ncol=3, nrow=1)
	ggsave(
		filename=paste0(plot,"combined_",gene,".",extension),
		device=extension,
		plot=arrangment,
		width=14,
		height=7
	)

	extension="png"
	arrangment = ggarrange(plots=plots_fixed,widths=c(4,4,4),ncol=3, nrow=1)
	ggsave(
		filename=paste0(plot,"combined_",gene,".",extension),
		device=extension,
		plot=arrangment,
		width=14,
		height=7
	)
}