source('~/necessary_packages.R')


skin_colors <- c(
  "Bulge cell" = "#E69F00",
  "Sebaceous gland" = "#56B4E9",
  "B keratinocyte" = "#009E73",
  "BS keratinocyte" = "#AD7700",
  "S keratinocyte" = "#D55E00",
  "G keratinocyte" = "#1C91D4",
  "Adipocyte" = "#007756",
  "Ccl21a+ Endothelial cell" = "#D5C711",
  "Cd300lg+ Endothelial cell" = "#005685",
  "Fibroblast" = "#A04700",
  "Langerhans cell" = "#B14380",
  "Melanocyte" = "#4D4D4D",
  "Schwan cell" = "#FFBE2D",
#  "Th2" = "#80C7EF",
  'Treg' = "#80C7EF",
  "Th17" = "#00F6B3",
  "Proliferating Th17" = "#F4EB71",
  "Th22" = "#06A5FF",
  "NK cell" = "#FF8320"
)


obj <- readRDS("/home2/sbjs0428/R/atopyNK/workspace/yesbatch/skin_ref2018_obj.RDS")



### NSK vs ATSK / ATSKTX vs ATSK volcano plot
deg <- FindMarkers(obj, group.by='Type', ident.1='ATSKTX', ident.2='ATSK') 
plot <- volcano(deg, title_ = 'DEGs between ATSKTX and ATSK')
ggsave('/home2/sbjs0428/R/atopyNK/plots/batch/skin/volcano/skin_ATSKTX_ATSK_volcano_8x10.png', plot, width=8, height=10, dpi=1200, units='in')

deg <- FindMarkers(obj, group.by='Type', ident.1='NSK', ident.2	='ATSK')
plot <- volcano(deg, title_ = 'DEGs between NSK and ATSK')
ggsave('/home2/sbjs0428/R/atopyNK/plots/batch/skin/volcano/skin_NSK_ATSK_volcano_8x10.png', plot, width=8, height=10, dpi=1200, units='in')



###cell type markers ###
genes <- c('Lgr5','Cd34','Nfatc1','Elovl5', 'Acsbg1', 'Far2', 'Acot5', 'Pnpla5', 'Krt14', 'Krt5', 'Krt15', 'Krt1', 'Krt10', 'Krt6a','Krt16', 'Krt79', 'S100a8', 'S100a9','Flg', 'Lor', 'Cfd', 'Car3', 'Cidec', 'Adig', 'Ptprb', 'Flt1', 'Vwf', 'Podxl','Ccl21a', 'Cd300lg', 'Dcn', 'Col3a1', 'Col1a2', 'Col1a1', 'Col6a3', 'Cd86', 'H2-M2', 'Cd207', 'Ltc4s', 'Dct', 'Tyrp1', 'Kit', 'Ptgds', 'Pax3', 'Syngr1', 'Mpz', 'Gatm', 'Mbp', 'Mal','Cd3e',  'Foxp3', 'Ctla4','Il2ra', 'Il7r', 'Tgfb1','Stat5a','Il10',  'Stat3','Ccr4', 'Ccr6','Il21','Il26','Rorc','Il17a', 'Il22','Tnf', 'Il13','Ahr','Bnc2',' Foxo4', 'Ccr10',  'Gzmb', 'Klrd1', 'Ptprcap', 'Il12rb', 'Prf1', 'Klrb1c', 'Klrk1')

genes <- genes[ genes %in% rownames(obj[['RNA']]$counts)]
png('/home2/sbjs0428/R/atopyNK/plots/batch/skin/heatmap/celltype_markers_heatmap_12x9.png', width=12, height=9, res=1200, units='in')
dittoHeatmap(obj, genes, annot.by='celltype', scaled.to.max=F, cluster_rows=F, annot.colors=skin_colors, heatmap.colors=colorRampPalette(c('white', 'red'))(35), slot='counts', scale='none', breaks=c(0:35), complex=F)
dev.off()




kera_markers <- c('Lgr5','Cd34','Nfatc1','Elovl5', 'Acsbg1', 'Far2', 'Acot5', 'Pnpla5', 'Krt14', 'Krt5', 'Krt15', 'Krt1', 'Krt10', 'Krt6a','Krt16', 'Krt79', 'S100a8', 'S100a9','Flg', 'Lor')

kera <- subset(obj, celltype %in% levels(obj$celltype)[1:6])
kera$celltype <- factor(kera$celltype, levels=c(as.vector(levels(obj$celltype)[1:6])))

png('/home2/sbjs0428/R/atopyNK/plots/batch/skin/skin_kera_dimplot_12x10.png', width=12, height=10, res=1200, units='in')
dittoDimPlot(kera, 'celltype', color.panel = unname(skin_colors)[1:6])
dev.off()


png('/home2/sbjs0428/R/atopyNK/plots/batch/skin/skin_dimplot_12x10.png', width=12, height=10, res=1200, units='in')
dittoDimPlot(obj, 'celltype', color.panel = unname(skin_colors))
dev.off()





## dittoHeatmap(kera, kera_markers, annot.by='celltype', scaled.to.max=F, cluster_rows=F, annot.colors=skin_colors, heatmap.colors=colorRampPalette(c('white', 'red'))(200), slot='counts', scale='none', breaks=c(0:200), complex=F )  -> can draw heatmap with raw counts, if complex set to F, breaks should be low or same to colors value 

#/home2/sbjs0428/R/atopyNK/plots/batch/skin/
png('/home2/sbjs0428/R/atopyNK/plots/batch/skin/heatmap/keratinocyte_markers_heatmap_12x9.png', width=12, height=9, res=1200, units='in')
dittoHeatmap(kera, kera_markers, annot.by='celltype', scaled.to.max=F, cluster_rows=F, annot.colors=skin_colors, heatmap.colors=colorRampPalette(c('white', 'red'))(50), slot='counts', scale='none', breaks=c(0:50), complex=F)
dev.off()

png('/home2/sbjs0428/R/atopyNK/plots/batch/skin/skin_kera_dotplot_8x10.png', width=8, height=10, res=1200, units='in')
dittoDotPlot(kera,  kera_markers, group.by='celltype')
dev.off()



png('/home2/sbjs0428/R/atopyNK/plots/batch/skin/skin_barplot_8x10.png', width=8, height=10, res=1200, units='in')
dittoBarPlot(obj, 'celltype', group.by='Type', retain.factor.levels=T, color.panel = unname(skin_colors))
dev.off()

source('/home2/sbjs0428/R/atopyNK/script/batch/subsampling.R')
obj_sampled <- subsampling()
png('/home2/sbjs0428/R/atopyNK/plots/batch/skin/skin_barplot_bulgenormalized_8x10.png', width=8, height=10, res=1200, units='in')
dittoBarPlot(obj_sampled, 'celltype', group.by='Type', retain.factor.levels = T, color.panel = unname(skin_colors))
dev.off()

png('/home2/sbjs0428/R/atopyNK/plots/batch/skin/skin_barplot_bulgenormalized_count_8x10.png', width=8, height=10, res=1200, units='in')
dittoBarPlot(obj_sampled, 'celltype', group.by='Type', retain.factor.levels = T, color.panel = unname(skin_colors), scale='count')
dev.off()


plist <- list()
plist[[1]] <- dittoDimPlot(obj, 'celltype', main='All samples', legend.show=F, color.panel = skin_colors)
plist[[2]] <- dittoDimPlot(obj, 'celltype', cells.use=obj$Type == 'NSK', main='NSK sample', legend.show=F, color.panel=skin_colors)
plist[[3]] <- dittoDimPlot(obj, 'celltype', cells.use=obj$Type == 'ATSK', main='ATSK sample', legend.show=F,color.panel=skin_colors)
plist[[4]] <- dittoDimPlot(obj, 'celltype', cells.use=obj$Type == 'ATSKTX', main='ATSKTX sample', legend.show=F, color.panel=skin_colors)
all_plots <- wrap_plots(plist, ncol=2, nrow=2)
ggsave('/home2/sbjs0428/R/atopyNK/plots/batch/skin/skin_dimplot_all-nsk-atsk-atsktx_10x10.png',all_plots, height=10, width=10, dpi=1200, units='in')

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
        cds <- learn_graph(cds, use_partition = F)
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

pseudo_boxplot <- function(cds, celltype_colors = NULL) {
data.pseudo <- as.data.frame(colData(cds))
p <- ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(celltype, monocle3_pseudotime, median), fill = celltype)) +
geom_boxplot() 
if (!is.null(celltype_colors)) {
	p <- p + scale_fill_manual(values = celltype_colors)
}
return(p)
}

### by type
pseudo_boxplot <- function(cds) {
data.pseudo <- as.data.frame(colData(cds))
ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(Type, monocle3_pseudotime, median), fill = Type)) +
geom_boxplot()
}




### version with cell number added
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






cds <- trajectory(kera, rootcell= c(10,11))
plot <- pseudo_boxplot(cds, skin_colors)
ggsave('/home2/sbjs0428/R/atopyNK/plots/batch/skin/skin_kera_pseudotime_boxplot_8x10.png', plot, height=10, width=8, dpi=1200, units='in')

plot <- pseudo_dimplot(cds)
ggsave('/home2/sbjs0428/R/atopyNK/plots/batch/skin/skin_kera_pseudo_10x10.png', plot, height=10, width=10, dpi=1200, units='in')

kera$pseudotime <- pseudotime(cds)
plot <- FeaturePlot(kera, 'pseudotime')
ggsave('/home2/sbjs0428/R/atopyNK/plots/batch/skin/skin_kera_pseudo_featureplot_10x10.png', plot, height=10, width=10, dpi=1200, units='in')


### for new version of boxplot)

temp_cds <- trajectory(kera, c(10,11))
plot <- pseudo_boxplot(temp_cds, celltype_colors = skin_colors, plot_title = 'Whole keratinocytes pseudotime boxplot')
ggsave('/home2/sbjs0428/R/atopyNK/plots/batch/skin/skin_kera_revised_pseudotime_boxplot_8x10.png', plot, height=10, width=8, dpi=1200, units='in')

for (type in unique(kera$Type)) {
        temp_obj <- subset(kera, Type == type)
        temp_cds <- trajectory(temp_obj, c(10,11))
        plot <- pseudo_boxplot(temp_cds, celltype_colors = skin_colors, plot_title = paste0(type, ' pseudotime boxplot'))
        ggsave(paste0('/home2/sbjs0428/R/atopyNK/plots/batch/skin/skin_kera_', type, '_revised_pseudotime_boxplot_8x10.png'), plot, height=10, width=8, dpi=1200, units='in')
}



####### original version ##########

for (type in unique(kera$Type)) {
        temp_obj <- subset(kera, Type == type)
        temp_cds <- trajectory(temp_obj, c(10,11))
	plot <- pseudo_boxplot(temp_cds, skin_colors)
	ggsave(paste0('/home2/sbjs0428/R/atopyNK/plots/batch/skin/skin_kera_', type, '_pseudotime_boxplot_8x10.png'), plot, height=10, width=8, dpi=1200, units='in')
#        plot <- pseudo_dimplot(temp_cds)
#        ggsave(paste0('/home2/sbjs0428/R/atopyNK/plots/batch/skin/skin_kera_', type, '_pseudo_10x10.png'), plot, height=10, width=10, dpi=1200, units='in')
#	temp_obj$pseudotime <- pseudotime(temp_cds)
#	plot <- FeaturePlot(temp_obj, 'pseudotime')
#	ggsave(paste0('/home2/sbjs0428/R/atopyNK/plots/batch/skin/skin_kera_', type, '_pseudo_featureplot_10x10.png'), plot, height=10, width=10, dpi=1200, units='in')
}




##markers between atsktx and atsk per celltype


##ontology patterns <- "leuko|cytokine|chemokine|wound|immune|skin|epi|interferon|macro|T cell"
patterns <- "inflam|leuko|cytokine|chemokine|wound|immune|skin|epithe|interferon|macrophage|T cell|kera|NK|interleukin|MHC|lympho|epidermal|natural killer|neutrophil|stimulus|potassium|chemotaxis|migration|infiltration"




############# Adipocyte ######################
#/home2/sbjs0428/R/atopyNK/ontology_output/batch/skin


adipo <- subset(obj, celltype == 'Adipocyte')
adipo_deg <- FindMarkers(adipo, group.by='Type', ident.1='ATSKTX', ident.2='ATSK')
adipo$celltype <- factor(adipo$celltype, levels = 'Adipocyte')

#plot <- volcano(adipo_deg, 'Adipocyte DEGs between ATSKTX and ATSK')#, ylim=c(0,120))
#ggsave('/home2/sbjs0428/R/atopyNK/plots/batch/skin/volcano/skin_adipo_volcano_8x10.png', plot, width=8, height=10, dpi=1200, units='in')


plot <- volcano(adipo_deg, 'Adipocyte DEGs between ATSKTX and ATSK', ylim=c(0,130))
ggsave('/home2/sbjs0428/R/atopyNK/plots/batch/skin/volcano/skin_adipo_volcano_8x10.png', plot, width=8, height=10, dpi=1200, units='in')


adipo_deg_filtered <- adipo_deg %>% filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.6)

if (nrow(adipo_deg_filtered) > 100) {
	adipo_deg_100 <- adipo_deg_filtered %>% arrange(desc(abs(pct.1-pct.2))) %>% slice_head(n=100)
	plot <- dittoHeatmap(adipo, rownames(adipo_deg_100), annot.by=c('celltype', 'Type'), order.by='Type', scaled.to.max=T, annot.colors=c("#007756", '#00BE67', "#377EB8", "#E41A1C"))
} else {
	plot <- dittoHeatmap(adipo, rownames(adipo_deg_filtered), annot.by=c('celltype', 'Type'), order.by='Type', scaled.to.max=T, annot.colors=c("#007756", '#00BE67', "#377EB8", "#E41A1C"))
}
ggsave('/home2/sbjs0428/R/atopyNK/plots/batch/skin/heatmap/adipo_markers_heatmap_10x13.png',plot, width=10, height=13, dpi=1200, units='in') 

#adipo_gse <- ontology(adipo_deg_atsktx_atsk, organism='mouse', num=15, title='Adipocyte DEGs', category.size=12, data.out=T)@result

adipo_gse_filtered <- ontology(adipo_deg %>% filter(p_val_adj < 0.05), organism='mouse', num=15, data.out=T, patterns=patterns)
write.xlsx(adipo_gse_filtered@result, '/home2/sbjs0428/R/atopyNK/ontology_output/batch/skin/adipo_ontology_output.xlsx', asTable=T)

#matching_descriptions <- grep(patterns, adipo_gse$Description, value=T)
#adipo_gse_filtered <- adipo_gse %>% filter(Description %in% matching_descriptions)

adipo_genes <- ''
for (i in adipo_gse_filtered$core_enrichment) {
	adipo_genes <- paste0(adipo_genes, i, sep='/')
}
adipo_genes <- unique(strsplit(adipo_genes, '/')[[1]])


plot <- only_plot_ontology(adipo_gse_filtered, num=30, title='Adipocyte DEGs filtered', category.size=10)
ggsave('/home2/sbjs0428/R/atopyNK/plots/batch/skin/skin_adipo_ontology_filtered_10x13.png',plot, width=10, height=13, dpi=1200, units='in')

plot <- dittoHeatmap(adipo, adipo_genes, annot.by=c('celltype', 'Type'), order.by='Type', scaled.to.max = T, annot.colors=c("#007756", '#00BE67', "#377EB8", "#E41A1C"))
ggsave('/home2/sbjs0428/R/atopyNK/plots/batch/skin/heatmap/adipo_markers_ontologyfiltered_heatmap_10x13.png',plot, width=10, height=13, dpi=1200, units='in')


########### Fibroblast ################


fibro <- subset(obj, celltype =='Fibroblast')
fibro$celltype <- factor(fibro$celltype, levels= 'Fibroblast')
fibro_deg <- FindMarkers(fibro, group.by='Type', ident.1='ATSKTX', ident.2='ATSK')
plot <- volcano(fibro_deg, 'Fibroblast DEGs between ATSKTX and ATSK', ylim=c(0,130))
ggsave('/home2/sbjs0428/R/atopyNK/plots/batch/skin/volcano/skin_fibro_volcano_8x10.png', plot, width=8, height=10, dpi=1200, units='in')

fibro_deg_filtered <- fibro_deg %>% filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.6) # %>% arrange(desc(abs(avg_log2FC))) %>% slice_head(n=100)
if (nrow(fibro_deg_filtered) > 100) {
        fibro_deg_100 <- fibro_deg_filtered %>% arrange(desc(abs(pct.1-pct.2))) %>% slice_head(n=100)
        plot <- dittoHeatmap(fibro, rownames(fibro_deg_100), annot.by=c('celltype', 'Type'), order.by='Type', scaled.to.max=T, annot.colors=c("#A04700", '#00BE67', "#377EB8", "#E41A1C"))
} else {
        plot <- dittoHeatmap(fibro, rownames(fibro_deg_filtered), annot.by=c('celltype', 'Type'), order.by='Type', scaled.to.max=T, annot.colors=c("#A04700", '#00BE67', "#377EB8", "#E41A1C"))
}
ggsave('/home2/sbjs0428/R/atopyNK/plots/batch/skin/heatmap/fibro_markers_heatmap_10x13.png',plot, width=10, height=13, dpi=1200, units='in')


fibro_gse_filtered <- ontology(fibro_deg_filtered, organism='mouse', num=12, title='Fibroblast DEGs', category.size=12, data.out=T, patterns = patterns)
write.xlsx(fibro_gse_filtered@result, '/home2/sbjs0428/R/atopyNK/ontology_output/batch/skin/fibro_ontology_output.xlsx', asTable=T)



fibro_genes <- ''
for (i in fibro_gse_filtered$core_enrichment) {
        fibro_genes <- paste0(fibro_genes, i, sep='/')
}
fibro_genes <- unique(strsplit(fibro_genes, '/')[[1]])

plot <- only_plot_ontology(fibro_gse_filtered, num=20, title='Fibroblast DEGs filtered', category.size=10)
ggsave('/home2/sbjs0428/R/atopyNK/plots/batch/skin/skin_fibro_ontology_filtered_10x13.png', plot, width=10, height=13, dpi=1200, units='in')

plot <- dittoHeatmap(fibro, fibro_genes, annot.by=c('celltype', 'Type'), order.by='Type', scaled.to.max = T, annot.colors=c("#A04700", '#00BE67', "#377EB8", "#E41A1C"))
ggsave('/home2/sbjs0428/R/atopyNK/plots/batch/skin/heatmap/fibro_markers_ontologyfiltered_heatmap_10x13.png',plot, width=10, height=13, dpi=1200, units='in')


############## Keratinocyte ##################


kera <- subset(obj, celltype %in% levels(obj$celltype)[3:6])
kera$celltype <- 'Keratinocyte'
kera$celltype <- factor(kera$celltype, levels ='Keratinocyte')
kera_deg <- FindMarkers(kera, group.by='Type', ident.1='ATSKTX', ident.2='ATSK')
plot <- volcano(kera_deg, 'Keratinocyte DEGs between ATSKTX and ATSK', ylim=c(0,130))
ggsave('/home2/sbjs0428/R/atopyNK/plots/batch/skin/volcano/skin_kera_volcano_8x10.png', plot, width=8, height=10, dpi=1200, units='in')

kera_deg_filtered <- kera_deg %>% filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.6) # %>% arrange(p_val_adj) %>% slice_head(n=100)
plot <- dittoHeatmap(kera, rownames(kera_deg_filtered), annot.by=c('celltype','Type'), order.by='Type', scaled.to.max=T, annot.colors=c("#D55E00", '#00BE67', "#377EB8", "#E41A1C"))
ggsave('/home2/sbjs0428/R/atopyNK/plots/batch/skin/heatmap/kera_markers_heatmap_10x13.png',plot, width=10, height=13, dpi=1200, units='in')

kera_gse <- ontology(kera_deg, organism='mouse', num=20, title='Keratinocyte DEGs', category.size=12, data.out=T)@result

kera_gse_filtered <- ontology(kera_deg_filtered, organism='mouse', num=12, title='Keratinocyte DEGs filtered', category.size=12, data.out=T, patterns = patterns)
write.xlsx(kera_gse_filtered@result, '/home2/sbjs0428/R/atopyNK/ontology_output/batch/skin/kera_ontology_output.xlsx', asTable=T)

kera_genes <- ''
for (i in kera_gse_filtered$core_enrichment) {
        kera_genes <- paste0(kera_genes, i, sep='/')
}
kera_genes <- unique(strsplit(kera_genes, '/')[[1]])



plot <- only_plot_ontology(kera_gse_filtered, num=20, title='Keratinocyte DEGs filtered', category.size=10)
ggsave('/home2/sbjs0428/R/atopyNK/plots/batch/skin/skin_kera_ontology_filtered_10x13.png',plot, width=10, height=13, dpi=1200, units='in')


kera_deg_ont <- kera_deg_filtered[kera_genes, ]
kera_deg_ont <- kera_deg_ont %>% arrange(desc(abs(pct.1-pct.2))) %>% slice_head(n=100)
plot <- dittoHeatmap(kera, rownames(kera_deg_ont), annot.by=c('celltype', 'Type'), order.by='Type', scaled.to.max = T, annot.colors=c("#CC79A7", '#00BE67', "#377EB8", "#E41A1C"))
ggsave('/home2/sbjs0428/R/atopyNK/plots/batch/skin/heatmap/kera_markers_ontologyfiltered_heatmap_10x15.png',plot, width=10, height=15, dpi=1200, units='in')


########### T cell subtypes #####################


tcells <- levels(obj$celltype)[c(14,15,17)]

for (t in tcells) {
	temp <- subset(obj, celltype == t)
	temp$celltype <- factor(temp$celltype, levels = t)
	deg_filtered <- FindMarkers(temp, group.by='Type', ident.1='ATSKTX', ident.2='ATSK') %>% filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.6)
	plot <- dittoHeatmap(temp, rownames(deg_filtered), annot.by=c('celltype', 'Type'), order.by='Type', scaled.to.max= T, annot.colors = c(skin_colors[[t]], '#00BE67', "#377EB8", "#E41A1C"))
	ggsave(paste0('/home2/sbjs0428/R/atopyNK/plots/batch/skin/heatmap/', t, '_markers_heatmap_10x13.png'), plot, widt=10, height=13, dpi=1200, units='in')
}



############## kera subtypes ################



keras <- levels(obj$celltype)[c(1,3,4,5,6)]

for (k in keras) {
	temp <- subset(obj, celltype == k)
	temp$celltype <- factor(temp$celltype, levels = k)
	deg_filtered <- FindMarkers(temp, group.by='Type', ident.1='ATSKTX', ident.2='ATSK') %>% filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.6)
	if (k %in% levels(obj$celltype)[c(1,3,5,6)]) {
		deg_filtered <- deg_filtered %>% arrange(desc(abs(pct.1-pct.2))) %>% slice_head(n=100)
	}
	plot <- dittoHeatmap(temp, rownames(deg_filtered), annot.by=c('celltype', 'Type'), order.by='Type', scaled.to.max= T, annot.colors = c(skin_colors[[k]], '#00BE67', "#377EB8", "#E41A1C"))
	ggsave(paste0('/home2/sbjs0428/R/atopyNK/plots/batch/skin/heatmap/', k, '_markers_heatmap_10x13.png'), plot, widt=10, height=13, dpi=1200, units='in')
}



















############# th2 or treg heatmap

plot <- dittoHeatmap(tcell, c('Ctla4', 'Foxp3', 'Il10', 'Stat5a', 'Stat6', 'Il4', 'Il13', 'Gata3'), annot.by=c('celltype', 'Type'),order.by='celltype',  scaled.to.max=F, annot.colors=c(unname(skin_colors[levels(tcell$celltype)]), '#00BE67', "#377EB8", "#E41A1C"), heatmap.colors=colorRampPalette(c('white', 'red'))(10), slot='counts', scale='none', breaks=c(0:10), complex=F, cluster_rows=F)




############# treg marker vlnplot

plot <- VlnPlot(tcell, c('Ctla4', 'Foxp3', 'Il10', 'Stat5a'), group.by='celltype', split.by='Type', add.noise=F, cols = c('#00BE67', "#377EB8", "#E41A1C"), ncol=2)



