#cell_markers <- list(
#  Naive_CD4_T_cell = c("Cd4", "Ptprc", "Sell", "Ccr7"),
#  Th1 = c("Tbx21", "Ifng", "Il12rb2", "Cxcr3", "Stat4"),
#  Th2 = c("Gata3", "Il4", "Il5", "Il13", "Il10"),
#  Treg = c("Foxp3", "Il2ra", "Ctla4", "Tnfrsf18", "Entpd1", "Nt5e"),
#  LTi = c("Rorc", "Cd4", "Il7r", "Thy1", "Lta", "Ltb", "Ccr6"),
#  CD4_CD8_Double_Negative_T_cell = c("Cd4", "Cd8a", "Cd8b1", "Trdc", "Trac"),
#  Th17 = c("Rorc", "Il17a", "Il17f", "Il21", "Il22"),
#  Th22 = c("Ahr", "Il22", "Il13", "Ccr6", "Ccr10"),
#  ILC = c("Id2", "Gata3", "Rorc", "Tbx21", "Il7r"),
#  gamma_delta_T_cell = c("Trdc", "Trdv1", 'ccr6', "Cd3e", "Cd3d")
#)
#
#
#    Follicular B_1     Follicular B_2       Marginal B_1       Marginal B_2 
#         "#E69F00"          "#56B4E9"          "#009E73"          "#F0E442" 
#          Plasma_1           Plasma_2   CD4- CD8- T cell   Naive CD4 T cell 
#         "#0072B2"          "#D55E00"          "#CC79A7"          "#666666" 
#               Th1                Th2               Treg gamma delta T cell 
#         "#AD7700"          "#1C91D4"          "#007756"          "#D5C711" 
#  Naive CD8 T cell      CD8 NK T cell            NK cell       Myeloid cell 
#         "#005685"          "#A04700"          "#B14380"          "#4D4D4D" 
#    Dendritic cell      Macrophage M1    Epithelial cell 
#         "#FFBE2D"          "#80C7EF"          "#00F6B3" 
#



lymph_colors <- c(
  "Follicular B_1" = "#E69F00",
  "Follicular B_2" = "#56B4E9",
  "Marginal B_1" = "#009E73",
  "Marginal B_2" = "#F0E442",
  "Plasma_1" = "#0072B2",
  "Plasma_2" = "#D55E00",
  "CD4- CD8- T cell" = "#CC79A7",
  "Naive CD4 T cell" = "#666666",
  "Th1" = "#AD7700",
  "Th2" = "#1C91D4",
  "Treg" = "#007756",
  "gamma delta T cell" = "#D5C711",
  "Naive CD8 T cell" = "#005685",
  "CD8 NK T cell" = "#A04700",
  "NK cell" = "#B14380",
  "Myeloid cell" = "#4D4D4D",
  "Dendritic cell" = "#FFBE2D",
  "Macrophage M1" = "#80C7EF",
  "Epithelial cell" = "#00F6B3"
)


source('~/necessary_packages.R')

obj <- readRDS('/home2/sbjs0428/R/atopyNK/workspace/nobatch/lymph_obj.RDS')


plot <- dittoDimPlot(obj, 'celltype') 
ggsave('/home2/sbjs0428/R/atopyNK/plots/nobatch/lymph/lymphnode_dimplot_12x10.png',plot, width=12, height=10, dpi=1200, units='in')

tcell <- subset(obj, celltype %in% levels(obj$celltype)[7:14])
tcell$celltype <- factor(tcell$celltype, levels=levels(obj$celltype)[7:14])



plist <- list()
plist[[1]] <- dittoDimPlot(obj, 'celltype', main='All samples', legend.show=F, color.panel=lymph_colors)
plist[[2]] <- dittoDimPlot(obj, 'celltype', cells.use=obj$Type == 'NLN', main='NLN sample', legend.show=F, color.panel=lymph_colors)
plist[[3]] <- dittoDimPlot(obj, 'celltype', cells.use=obj$Type == 'ATLN', main='ATLN sample', legend.show=F, color.panel=lymph_colors)
plist[[4]] <- dittoDimPlot(obj, 'celltype', cells.use=obj$Type == 'ATLNTX', main='ATLNTX sample', legend.show=F, color.panel=lymph_colors)
all_plots <- wrap_plots(plist, ncol=2, nrow=2)
ggsave('/home2/sbjs0428/R/atopyNK/plots/nobatch/lymph/lymph_dimplot_all-nln-atln-atlntx_10x10.png',all_plots, height=10, width=10, dpi=1200, units='in')

rm(plist, all_plots)

plot <- dittoBarPlot(obj, 'celltype', group.by='Type', retain.factor.levels=T)
ggsave('/home2/sbjs0428/R/atopyNK/plots/nobatch/lymph/lymphnode_barplot_8x10.png', plot,width=8, height=10, dpi=1200, units='in')

rm(plot)

tcell <- subset(obj, celltype %in% levels(obj$celltype)[7:14])
cd4 <- subset(obj, celltype %in% levels(obj$celltype)[7:12])
cd8 <- subset(obj, celltype %in% c('CD4- CD8- T cell', levels(obj$celltype)[13:14]))

#tcell_atln <- subset(tcell, Type == 'ATLN')
#tcell_nln <- subset(tcell, Type == 'NLN')
#tcell_atlntx <- subset(tcell, Type =='ATLNTX')


#Th2 deg heatmap and dotplot
th2_deg <- FindMarkers(obj, group.by='celltype', ident.1='Th2', only.pos=T) %>% filter(p_val_adj < 0.05) %>% filter(avg_log2FC > 5)

th2_markers <- union(c("Gata3", "Il4", "Il5", "Il13", "Il10"), rownames(th2_deg))

plot <- dittoHeatmap(obj, th2_markers, annot.by='Type', scaled.to.max = T)
ggsave('/home2/sbjs0428/R/atopyNK/plots/nobatch/lymph/heatmap/lymph_th2_degfromWhole_heatmap_10x13.png',plot, width=10, height=13, units='in', dpi=1200)



trajectory <- function(obj, rootcell) {
        cds <- as.cell_data_set(obj)
        fData(cds)$gene_short_name <- rownames(fData(cds))
        reacreate.partition <- c(rep(1,length(cds@colData@rownames)))
        names(reacreate.partition) <- cds@colData@rownames
        reacreate.partition <- as.factor(reacreate.partition)
        cds@clusters$UMAP$partitions <- reacreate.partition
        list_cluster <- obj@active.ident
        cds@clusters$UMAP$clusters <- list_cluster
        cds@int_colData@listData$reducedDims$UMAP <- obj@reductions$umap@cell.embeddings
        cds <- learn_graph(cds, use_partition = FALSE)
#       plot_cells(cds,
#           color_cells_by = 'celltype',
#           label_groups_by_cluster = FALSE,
#           label_branch_points = FALSE,
#           label_roots = FALSE,
#           label_leaves = FALSE,
#           group_label_size = 5)

        cds <- order_cells(cds, reduction_method = 'UMAP', root_cells = colnames(cds)[ cds$seurat_clusters == rootcell])
        cds$monocle3_pseudotime <- pseudotime(cds)
        return(cds)
}

pseudo_dimplot <- function(cds) {
plot_cells(cds,
color_cells_by = 'pseudotime',
label_groups_by_cluster = FALSE,
label_branch_points = FALSE,
label_roots = FALSE,
label_leaves = FALSE)
}

pseudo_boxplot <- function(cds) {
data.pseudo <- as.data.frame(colData(cds))
ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(celltype, monocle3_pseudotime, median), fill = celltype)) +
geom_boxplot()
}


## version with cell number added
pseudo_boxplot <- function(cds, celltype_colors = NULL, plot_title = "Pseudotime Boxplot") {
  data.pseudo <- as.data.frame(colData(cds))

  p <- ggplot(data.pseudo, aes(x = monocle3_pseudotime, y = reorder(celltype, monocle3_pseudotime, median), fill = celltype)) +
    geom_violin(trim = FALSE, scale = "width", alpha = 0.6) +  # Adds a violin plot to show the distribution
    geom_boxplot(width = 0.2, outlier.shape = NA) +  # Adds a boxplot inside the violin plot
    stat_summary(fun.data = function(y) data.frame(y = max(y), label = length(y)),
                 geom = "text", size = 3, hjust = -0.3) +  # Adds the number of cells as text
    labs(x = "Pseudotime", y = "Cell Type", title = plot_title) +  # Adds title to the plot
    theme_minimal() +
    theme(
      legend.position = "none",  # Hides the legend for fill
      panel.background = element_rect(fill = "white", color = NA),  # Set panel background color to white
      plot.background = element_rect(fill = "white", color = NA),  # Set plot background color to white
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank()   # Remove minor grid lines
    )

  if (!is.null(celltype_colors)) {
    p <- p + scale_fill_manual(values = celltype_colors)  # Use the provided colors for cell types
  }

  print(p)
}







#################### T cell pseudotime

cds_cd4t <- trajectory(cd4, 12)
cd4$pseudotime <- pseudotime(cds_cd4t)
plot <- pseudo_dimplot(cds_cd4t)
ggsave('/home2/sbjs0428/R/atopyNK/plots/nobatch/lymph/pseudo_cd4/lymph_cd4t_pseudo_dimplot_10x10.png', plot, width=10, height=10, dpi=1200, units='in')
plot <- FeaturePlot(cd4, 'pseudotime')
ggsave('/home2/sbjs0428/R/atopyNK/plots/nobatch/lymph/pseudo_cd4/lymph_cd4t_pseudo_featureplot_10x10.png', plot, width=10, height=10, dpi=1200, units='in')



### for revised boxplot version
cd4$celltype <- factor(cd4$celltype, levels = levels(obj$celltype)[7:12])
temp_cds <- trajectory(cd4, 12)
plot <- pseudo_boxplot(temp_cds, lymph_colors, 'Whole CD4 T cell pseudotime boxplot')
ggsave('/home2/sbjs0428/R/atopyNK/plots/nobatch/lymph/pseudo_cd4/lymph_cd4_revised_pseudotime_boxplot_8x10.png', width=8, height=10, dpi=1200, units='in')


for (type in unique(cd4$Type)) {
    temp_cd4 <- subset(cd4, Type == type)
    temp_cds <- trajectory(temp_cd4, 12)
    plot <- pseudo_boxplot(temp_cds, lymph_colors, paste0(type, ' CD4 T cell pseudotime boxplot'))
    ggsave(paste0('/home2/sbjs0428/R/atopyNK/plots/nobatch/lymph/pseudo_cd4/lymph_cd4_', type, '_revised_pseudotime_boxplot_8x10.png'), width=8, height=10, dpi=1200, units='in')
} 



tcell <- subset(obj, celltype %in% levels(obj$celltype)[7:14])
tcell$celltype <- factor(tcell$celltype, levels=levels(obj$celltype)[7:14])

temp_cds <- trajectory(tcell, 12)
plot <- pseudo_boxplot(temp_cds, lymph_colors, 'Whole T cell pseudotime boxplot')
ggsave('/home2/sbjs0428/R/atopyNK/plots/nobatch/lymph/pseudo_tcell/lymph_tcell_revised_pseudotime_boxplot_8x10.png', width=8, height=10, dpi=1200, units='in')


for (type in unique(tcell$Type)) {
    temp_tcell <- subset(tcell, Type == type)
    temp_cds <- trajectory(temp_tcell, 12)
    plot <- pseudo_boxplot(temp_cds, lymph_colors, paste0(type, ' T cell pseudotime boxplot'))
    ggsave(paste0('/home2/sbjs0428/R/atopyNK/plots/nobatch/lymph/pseudo_tcell/lymph_tcell_', type, '_revised_pseudotime_boxplot_8x10.png'), width=8, height=10, dpi=1200, units='in')
} 









for (type in unique(cd4$Type)) {
    temp_cd4 <- subset(cd4, Type == type)
    temp_cds <- trajectory(temp_cd4, 12)
    temp_cd4$pseudotime <- pseudotime(temp_cds)
    
    # Plot the pseudotime dimensional plot
    plot <- pseudo_dimplot(temp_cds)
    ggsave(paste0('/home2/sbjs0428/R/atopyNK/plots/nobatch/lymph/pseudo_cd4/lymph_cd4_', type, '_pseudo_dimplot_10x10.png'), plot, width=10, height=10, dpi=1200, units='in')
    
    # Plot the FeaturePlot for pseudotime
    plot <- FeaturePlot(temp_cd4, 'pseudotime')
    ggsave(paste0('/home2/sbjs0428/R/atopyNK/plots/nobatch/lymph/pseudo_cd4/lymph_cd4_', type, '_pseudo_featureplot_10x10.png'), plot, width=10, height=10, dpi=1200, units='in')}



cds_tcell <- trajectory(tcell, 12)
tcell$pseudotime <- pseudotime(cds_tcell)
plot <- pseudo_dimplot(cds_tcell)
ggsave('/home2/sbjs0428/R/atopyNK/plots/nobatch/lymph/lymph_tcell_pseudo_dimplot_10x10.png', plot, width=10, height=10, dpi=1200, units='in')
plot <- FeaturePlot(tcell, 'pseudotime')
ggsave('/home2/sbjs0428/R/atopyNK/plots/nobatch/lymph/lymph_tcell_pseudo_featureplot_10x10.png', plot, width=10, height=10, dpi=1200, units='in')



for (type in unique(tcell$Type)) {
	temp_obj <- subset(tcell, Type == type)
	temp_obj$celltype <- factor(temp_obj$celltype, levels = levels(obj$celltype)[7:14])
	temp_cds <- trajectory(temp_obj, 12)
	plot <- pseudo_dimplot(temp_cds)
	ggsave(paste0('/home2/sbjs0428/R/atopyNK/plots/nobatch/lymph/pseudo_tcell/lymph_tcell_', type, '_pseudo_dimplot_10x10.png'), plot, width=10, height=10, dpi=1200, units='in')
	temp_obj$pseudotime <- pseudotime(temp_cds)
	plot <- FeaturePlot(temp_obj, 'pseudotime')
	ggsave(paste0('/home2/sbjs0428/R/atopyNK/plots/nobatch/lymph/pseudo_tcell/lymph_tcell_', type, '_pseudo_featureplot_10x10.png'), plot, width=10, height=10, dpi=1200, units='in')
}
	


#cd4$pseudotime <- pseudotime(cds)
#Idents(obj) <- obj$celltype
#return(obj)


tcells <- c('tcell', 'cd4', 'cd8')
samples <- c('all', 'nln', 'atln', 'atlntx')


for (name in tcells) {
	for (sample in samples) {
		if (sample == 'all') {
			type <- get(name)
		} else {
			type <- subset(get(name), subset = Type == toupper(sample))
		}
		cds <- trajectory(type, rootcell=12)
		plot <- pseudo_boxplot(cds)
		dir_path <- paste0('/home2/sbjs0428/R/atopyNK/plots/nobatch/lymph/', 'pseudo_', name)
    		if (!dir.exists(dir_path)) {
      			dir.create(dir_path, recursive = TRUE)
    		}
		if ( sample == 'all') {
			ggsave(paste0('/home2/sbjs0428/R/atopyNK/plots/nobatch/lymph/', 'pseudo_', name, '/', 'lymph_',name,'_pseudo_boxplot_8x10.png'), plot, width=8, height=10, dpi=1200, units='in')
		} else {
			ggsave(paste0('/home2/sbjs0428/R/atopyNK/plots/nobatch/lymph/', 'pseudo_', name, '/', 'lymph_',name, '_',sample, '_pseudo_boxplot_8x10.png'), plot, width=8, height=10, dpi=1200, units='in')
		}
	}
}

########################################################################

#annot.colors=c("#007756", '#00BE67', "#377EB8", "#E41A1C"))



## deg heatmap by celltype

patterns <- "inflam|leuko|cytokine|chemokine|wound|immune|skin|epithe|interferon|macrophage|T cell|kera|NK|interleukin|MHC|lympho|epidermal|natural killer|neutrophil|stimulus"


## cd4 t cells

cd4t <- subset(obj, celltype %in% levels(obj$celltype)[8:12])
cd4t$celltype <- factor("CD4 T cell", levels = 'CD4 T cell')
cd4t_deg <- FindMarkers(cd4t, group.by='Type', ident.1='ATLNTX', ident.2='ATLN') %>% filter(p_val_adj < 0.05)
cd4t_deg_filtered <- cd4t_deg %>% filter(p_val_adj < 0.05, abs(avg_log2FC) > 1)
plot <- dittoHeatmap(cd4t, rownames(cd4t_deg_filtered), annot.by=c('celltype', 'Type'), order.by='Type', scaled.to.max=T, annot.colors=c("#666666", '#00BE67', "#377EB8", "#E41A1C"))
ggsave('/home2/sbjs0428/R/atopyNK/plots/nobatch/lymph/heatmap/cd4t_markers_heatmap_10x13.png',plot, width=10, height=13, dpi=1200, units='in')

cd4t_deg <- FindMarkers(cd4t, group.by='Type', ident.1='ATLNTX', ident.2='ATLN')
cd4t_gse <- ontology(cd4t_deg, organism='mouse', data.out=T)
plot <- only_plot_ontology(cd4t_gse, num=15, title='CD4 T cell DEGs', category.size=10)
ggsave('/home2/sbjs0428/R/atopyNK/plots/nobatch/lymph/lymph_cd4t_ontology_10x13.png', plot, width=10, height=13, dpi=1200, units='in')



cd4t_gse_filtered <- ontology(cd4t_deg, organism='mouse', data.out=T, patterns=patterns)
write.xlsx(cd4t_gse_filtered@result, '/home2/sbjs0428/R/atopyNK/ontology_output/nobatch/lymph/cd4t_ontology_output.xlsx', asTable=T)
cd4t_genes <- ''
for (i in cd4t_gse_filtered$core_enrichment) {
	cd4t_genes <- paste0(cd4t_genes, i, sep='/')
}
cd4t_genes <- unique(strsplit(cd4t_genes, '/')[[1]])
plot <- only_plot_ontology(cd4t_gse_filtered, num=30, title='CD4 T cell DEGs filtered', category.size=10)
ggsave('/home2/sbjs0428/R/atopyNK/plots/nobatch/lymph/lymph_cd4t_ontology_filtered_10x13.png',plot, width=10, height=13, dpi=1200, units='in')
plot <- dittoHeatmap(cd4t, cd4t_genes, annot.by=c('celltype', 'Type'), order.by='Type', scaled.to.max = T, annot.colors=c("#666666", '#00BE67', "#377EB8", "#E41A1C"))
ggsave('/home2/sbjs0428/R/atopyNK/plots/nobatch/lymph/heatmap/cd4t_markers_ontologyfiltered_heatmap_10x13.png',plot, width=10, height=13, dpi=1200, units='in')




## plasma 

plasma <- subset(obj, celltype %in% levels(obj$celltype)[5:6])
plasma$celltype <- factor('Plasma cell', levels = 'Plasma cell')
plasma_deg <- FindMarkers(plasma, group.by='Type', ident.1='ATLNTX', ident.2='ATLN') %>% filter(p_val_adj < 0.05)
plasma_deg_filtered <- plasma_deg %>% filter(p_val_adj < 0.05, abs(avg_log2FC) > 1)
plot <- dittoHeatmap(plasma, rownames(plasma_deg_filtered), annot.by=c('celltype', 'Type'), order.by='Type', scaled.to.max=T, annot.colors=c("#D55E00", '#00BE67', "#377EB8", "#E41A1C"))
ggsave('/home2/sbjs0428/R/atopyNK/plots/nobatch/lymph/heatmap/plasma_markers_heatmap_10x13.png',plot, width=10, height=13, dpi=1200, units='in')

plasma_deg <- FindMarkers(plasma, group.by='Type', ident.1='ATLNTX', ident.2='ATLN')
plasma_gse <- ontology(plasma_deg, organism='mouse', data.out=T)
plot <- only_plot_ontology(plasma_gse, num=15, title='Plasma cell DEGs', category.size=10)
ggsave('/home2/sbjs0428/R/atopyNK/plots/nobatch/lymph/lymph_plasma_ontology_10x13.png', plot, width=10, height=13, dpi=1200, units='in')



plasma_gse_filtered <- ontology(plasma_deg, organism='mouse', data.out=T, patterns=patterns)
write.xlsx(plasma_gse_filtered@result, '/home2/sbjs0428/R/atopyNK/ontology_output/nobatch/lymph/plasma_ontology_output.xlsx', asTable=T)
plasma_genes <- ''
for (i in plasma_gse_filtered$core_enrichment) {
        plasma_genes <- paste0(plasma_genes, i, sep='/')
}
plasma_genes <- unique(strsplit(plasma_genes, '/')[[1]])
plot <- only_plot_ontology(plasma_gse_filtered, num=30, title='Plasma cell DEGs filtered', category.size=10)
ggsave('/home2/sbjs0428/R/atopyNK/plots/nobatch/lymph/lymph_plasma_ontology_filtered_10x13.png',plot, width=10, height=13, dpi=1200, units='in')
plot <- dittoHeatmap(plasma, plasma_genes, annot.by=c('celltype', 'Type'), order.by='Type', scaled.to.max = T, annot.colors=c("#D55E00", '#00BE67', "#377EB8", "#E41A1C"))
ggsave('/home2/sbjs0428/R/atopyNK/plots/nobatch/lymph/heatmap/plasma_markers_ontologyfiltered_heatmap_10x13.png',plot, width=10, height=13, dpi=1200, units='in')




## macro

macro <- subset(obj, celltype == 'Macrophage M1')
macro$celltype <- factor(macro$celltype, levels = unique(macro$celltype))
macro_deg <- FindMarkers(macro, group.by='Type', ident.1='ATLNTX', ident.2='ATLN') %>% filter(p_val_adj < 0.05)
macro_deg_filtered <- macro_deg %>% filter(p_val_adj < 0.05, abs(avg_log2FC) > 1)
plot <- dittoHeatmap(macro, rownames(macro_deg_filtered), annot.by=c('celltype', 'Type'), order.by='Type', scaled.to.max=T, annot.colors=c("#80C7EF", '#00BE67', "#377EB8", "#E41A1C"))
ggsave('/home2/sbjs0428/R/atopyNK/plots/nobatch/lymph/heatmap/macro_markers_heatmap_10x13.png',plot, width=10, height=13, dpi=1200, units='in')


macro_deg <- FindMarkers(macro, group.by='Type', ident.1='ATLNTX', ident.2='ATLN')
macro_gse <- ontology(macro_deg, organism='mouse', data.out=T)
plot <- only_plot_ontology(macro_gse, num=15, title='Macrophage M1 DEGs', category.size=10)
ggsave('/home2/sbjs0428/R/atopyNK/plots/nobatch/lymph/lymph_macro_ontology_10x13.png', plot, width=10, height=13, dpi=1200, units='in')



macro_gse_filtered <- ontology(macro_deg, organism='mouse', data.out=T, patterns=patterns)
write.xlsx(macro_gse_filtered@result, '/home2/sbjs0428/R/atopyNK/ontology_output/nobatch/lymph/macro_ontology_output.xlsx', asTable=T)
macro_genes <- ''
for (i in macro_gse_filtered$core_enrichment) {
        macro_genes <- paste0(macro_genes, i, sep='/')
}
macro_genes <- unique(strsplit(macro_genes, '/')[[1]])
plot <- only_plot_ontology(macro_gse_filtered, num=30, title='Macrophage M1 DEGs filtered', category.size=10)
ggsave('/home2/sbjs0428/R/atopyNK/plots/nobatch/lymph/lymph_macro_ontology_filtered_10x13.png',plot, width=10, height=13, dpi=1200, units='in')
plot <- dittoHeatmap(macro, macro_genes, annot.by=c('celltype', 'Type'), order.by='Type', scaled.to.max = T, annot.colors=c("#80C7EF", '#00BE67', "#377EB8", "#E41A1C"))
ggsave('/home2/sbjs0428/R/atopyNK/plots/nobatch/lymph/heatmap/macro_markers_ontologyfiltered_heatmap_10x13.png',plot, width=10, height=13, dpi=1200, units='in')







## B cell

bcell <- subset(obj, celltype %in% levels(obj$celltype)[1:4])
bcell$celltype <- factor('B cell', levels = 'B cell')
bcell_deg <- FindMarkers(bcell, group.by='Type', ident.1='ATLNTX', ident.2='ATLN')
bcell_deg_filtered <- bcell_deg %>% filter(p_val_adj < 0.05, abs(avg_log2FC) > 1)


plot <- dittoHeatmap(bcell, rownames(bcell_deg_filtered), annot.by=c('celltype', 'Type'), order.by='Type', scaled.to.max=T, annot.colors=c("#E69F00", '#00BE67', "#377EB8", "#E41A1C"))
ggsave('/home2/sbjs0428/R/atopyNK/plots/nobatch/lymph/heatmap/bcell_markers_heatmap_10x13.png', plot, width=10, height=13, dpi=1200, units='in')



bcell_deg_filtered <- bcell_deg %>% filter(p_val_adj < 0.05) ## filter only p values
plot <- ontology(bcell_deg_filtered, organism='mouse', num=15, title='B cell DEGs', category.size=10)
