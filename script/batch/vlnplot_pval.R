library(patchwork)
library(Seurat)
library(ggplot2)
library(ggsignif)

#my_comparisons = list(c('ATSKTX', 'ATSK'), c('NSK','ATSK'))
##skin, avg_log2FC absolute value > 1 
#
#name <- list(1,2,3,4,5,6,7,8,9,10,11,12)
#names(name) <- c('Hrnr', 'Stfa2', 'Flg', 'Lor', 'Ivl', 'Csta1', 'Cnfn','Il18', 'Cysrt1', 'Serpinb12', 'Il15', 'Il13')
#
##name <- list(1,2,3,4,5,6,7,8,9)
##names(name) <- c('Stfa2', 'Flg', 'Lor', 'Ivl', 'Csta1', 'Cnfn', 'Cysrt1', 'Serpinb12', 'Il15')
#
#plot_list <- list()
#
#for (i in 1:12){
#	plot = VlnPlot(obj, names(name[i]), group.by='Type', add.noise=F, cols=c('#00BE67', "#377EB8", "#E41A1C")) + geom_signif(comparisons= my_comparisons, map_signif_level= function(p) sprintf("*p = %.2g", p), step_increase=0.15) + ylim(0, 9)
#	#assign(paste0('p', i) ,plot, envir=.GlobalEnv)
#	plot_list[[i]] <- plot
# }
#combined_plot <- wrap_plots(plot_list, ncol=4, nrow=3)
#ggsave('/home2/sbjs0428/R/atopyNK/plots/nobatch/skin/skin_ontology_atopyrelated_vlnplot_18x13.png', combined_plot, width = 18, height = 13, dpi = 1200, units = 'in')
##ggsave('plots/skin/skin_ontology_atopyrelated_vlnplot_15x13.png', combined_plot, width = 15, height = 13, dpi = 1200, units = 'in')


#lymph, avg_log2FC absolute value 0.5 adjusted 

# ontology terms related to atopy
num <- read.table('deg_atlntx_atln_ontology_output_related_num.txt')
name <- vector()
for ( i in num) {
	name <- c(name, i)
}

categories <- data.frame(gse@result)[name, 'Description']
#dotplot(gse, showCategory=categories, orderBy = "setSize", split=".sign") +  facet_grid(.~.sign)


gene_list <- vector()
for(index in num) {
g <- strsplit(gse@result[index, 'core_enrichment'], '/')[[1]]
gene_list <- unique(c(gene_list, g))
}


filtered_ont_output <-read.xlsx('ontology_output/lymph/Filtered_Immune_Related_Genes_lymphnode.xlsx')
gene_list <- vector()
for(index in 1:nrow(filtered_ont_output)) {
	genes <- strsplit(filtered_ont_output$Genes[index], ', ')[[1]]
	gene_list <- unique(c(gene_list, genes))
}

selected_genes <- vector()
for(line in readLines('lymph/atopy_related_degs.txt')) {
	g <- strsplit(line, ' ')[[1]][1]
	selected_genes <- c(selected_genes, g)
}

selected_genes <- c('Ifng', 'Gfi1', 'Il10', 'Slpi', 'Serpinb9', 'Prg2', 'Cxcr3', 'Ccr5', 'Fcer2a', 'Cxcl9', 'Ier3', 'Il12rb2', 'Irf4', 'Ccr2', 'Osm')

name <- list(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)
names(name) <- selected_genes

plot_list <- list()

for (i in 1:15){
        plot = VlnPlot(th2, names(name[i]), group.by='Type', add.noise=F, cols=c('#00BE67', "#377EB8", "#E41A1C")) + geom_signif(comparisons= my_comparisons, map_signif_level= function(p) sprintf("*p = %.2g", p), step_increase=0.15) + ylim(0, 9)
        #assign(paste0('p', i) ,plot, envir=.GlobalEnv)
        plot_list[[i]] <- plot
 }
combined_plot <- wrap_plots(plot_list, ncol=5, nrow=3)
ggsave('/home2/sbjs0428/R/atopyNK/plots/nobatch/lymph/lymph_ontology_atopyrelated_vlnplot_23x12.png', combined_plot, width = 23, height = 12, dpi = 1200, units = 'in')


VlnPlot(obj, g, group.by='Type', add.noise=F, cols=c('#00BE67', "#377EB8", "#E41A1C")) + geom_signif(comparisons= my_comparisons, map_signif_level= function(p) sprintf("*p = %.2g", p), step_increase=0.15) + ylim(0, 9)




name <- lapply(1:4, function(i) i)

names(name) <- c('Ctla4', 'Foxp3', 'Il10', 'Stat5a')

plot_list <- list()

for (i in 1:4) {
	plot = VlnPlot(tcell, names(name[i]), group.by='Type', add.noise=F, cols=c('#00BE67', "#377EB8", "#E41A1C")) + geom_signif(comparisons= my_comparisons, map_signif_level= function(p) sprintf("*p = %.2g", p), step_increase=0.15) + ylim(0, 9)
        #assign(paste0('p', i) ,plot, envir=.GlobalEnv)
        plot_list[[i]] <- plot
 }

combined_plot <- wrap_plots(plot_list, ncol=2, nrow=2)
