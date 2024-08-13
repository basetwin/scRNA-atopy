source('~/necessary_packages.R')

#obj <- readRDS("/home2/sbjs0428/R/atopyNK/workspace/yesbatch/skin_ref2018_obj.RDS")


#DimPlot drawing
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







Idents(obj) <- obj$celltype
objlist <- SplitObject(obj, split.by = 'Type')



for (name in names(objlist)) {
	object <- objlist[[name]]
	Idents(object) <- object$celltype
	labels <- Idents(object)
	meta <- data.frame(labels = labels, row.names = names(labels))
	cellchat <- createCellChat(object, group.by='labels', assay='RNA', meta = meta)
	cellchat <- setIdent(cellchat, ident.use = 'labels')
	cellchat@DB <- CellChatDB.mouse
	cellchat <- subsetData(cellchat)
	cellchat <- identifyOverExpressedGenes(cellchat)
	cellchat <- identifyOverExpressedInteractions(cellchat)
	cellchat <- projectData(cellchat, PPI.mouse)
	cellchat <- computeCommunProb(cellchat, type = "triMean")
	cellchat <- filterCommunication(cellchat, min.cells = 10)
	# save
	cellchat <- computeCommunProbPathway(cellchat)
	
	df.net <- subsetCommunication(cellchat)
	cellchat <- computeCommunProbPathway(cellchat)
	cellchat <- aggregateNet(cellchat)
	cellchat <- netAnalysis_computeCentrality(cellchat, slot.name='netP')
	saveRDS(cellchat, paste0('/home2/sbjs0428/R/atopyNK/cellchat/batch/skin/data/',name, '_cellchat.RDS'))
}


at <- readRDS('/home2/sbjs0428/R/atopyNK/cellchat/batch/skin/data/ATSK_cellchat.RDS')
atx <- readRDS('/home2/sbjs0428/R/atopyNK/cellchat/batch/skin/data/ATSKTX_cellchat.RDS')
object.list <- list('ATSK' = at, 'ATSKTX' = atx)
saveRDS(object.list, '/home2/sbjs0428/R/atopyNK/cellchat/batch/skin/data/objlist.RDS')
cellchat <- mergeCellChat(object.list, add.names = names(object.list))



#object <- obj

#Idents(object) <- object$celltype
#        labels <- Idents(object)
#        meta <- data.frame(labels = labels, row.names = names(labels))
#        cellchat <- createCellChat(object, group.by='labels', assay='RNA', meta = meta)
#        cellchat <- setIdent(cellchat, ident.use = 'labels')
#        cellchat@DB <- CellChatDB.mouse
#        cellchat <- subsetData(cellchat)
#        cellchat <- identifyOverExpressedGenes(cellchat)

#        cellchat <- identifyOverExpressedInteractions(cellchat)
#        cellchat <- projectData(cellchat, PPI.mouse)
#        cellchat <- computeCommunProb(cellchat, type = "triMean")
#        cellchat <- filterCommunication(cellchat, min.cells = 10)
#        # save
#        df.net <- subsetCommunication(cellchat)
#        cellchat <- computeCommunProbPathway(cellchat)
#        cellchat <- aggregateNet(cellchat)
#        saveRDS(cellchat, 'skin_cellchat.RDS')

#color
#ATSK -> "#377EB8"
#ATSKTX -> "#E41A1C"
#c("#377EB8", "#E41A1C")

#further step after merging atsk and atsktx 

#before run the code, python should be loaded
reticulate::use_python('~/miniconda3/envs/monocle3/bin/python', required=T)

#cellchat <- netAnalysis_computeCentrality(cellchat, slot.name='netP')

cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")

cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")



cellchat@idents$joint <- factor(cellchat@idents$joint, levels= levels(obj$celltype))

saveRDS(cellchat, '/home2/sbjs0428/R/atopyNK/cellchat/batch/skin/data/cellchat_merged_at_atx.RDS')

pathway.union <- union(object.list[[1]]@netP$pathways, object.list[[2]]@netP$pathways)



library(ComplexHeatmap)
png('/home2/sbjs0428/R/atopyNK/plots/batch/skin/cellchat/atsk_atsktx_pathway_barplot_8x18.png', width=8, height=18, units='in', res=1200)
rankNet(cellchat, mode='comparison', measure='weight', sources.use=NULL, targets.use=NULL, stacked=T, do.stat=T, color.use=c("#377EB8", "#E41A1C"))
dev.off()


source('/home2/sbjs0428/R/atopyNK/cellchat/batch/skin/script/netVisual_heatmap.R')

png('/home2/sbjs0428/R/atopyNK/plots/batch/skin/cellchat/atsk_atsktx_differential_weight_heatmap_12x10.png', width=12, height=10, units='in', res=1200)
ht <- netVisual_heatmap(cellchat, measure='weight', title.name='Differential interactions-weight ATSK vs ATSKTX', color.use=skin_colors)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right",
     padding = unit(c(5, 5, 2, 0), "mm"))  # top, right, bottom, left padding
dev.off()

png('/home2/sbjs0428/R/atopyNK/plots/batch/skin/cellchat/atsk_atsktx_differential_count_heatmap_12x10.png', width=12, height=10, units='in', res=1200)
ht <- netVisual_heatmap(cellchat, measure='count', title.name='Differential interactions-count ATSK vs ATSKTX', color.use=skin_colors)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right",
     padding = unit(c(5, 5, 2, 0), "mm"))  # top, right, bottom, left padding
dev.off()

png('/home2/sbjs0428/R/atopyNK/plots/batch/skin/cellchat/atsk_atsktx_outgoing_heatmap_11x16.png', width=11, height=16, units='in', res=1200)

p1 <- netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[1], width = 10, height = 33, color.use=unname(skin_colors[levels(object.list[[1]]@meta$labels)]))

p2 <- netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[2], width = 10, height = 33, color.use=unname(skin_colors[levels(object.list[[2]]@meta$labels)]))

draw(p1+p2, ht_gap=unit(0.5, 'cm'))
dev.off()

png('/home2/sbjs0428/R/atopyNK/plots/batch/skin/cellchat/atsk_atsktx_incoming_heatmap_11x16.png', width=11, height=16, units='in', res=1200)
p1 <- netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[1], width = 10, height = 33, color.use=unname(skin_colors[levels(object.list[[1]]@meta$labels)]))
p2 <- netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[2], width = 10, height = 33, color.use=unname(skin_colors[levels(object.list[[2]]@meta$labels)]))
draw(p1+p2, ht_gap=unit(0.5, 'cm'))
dev.off()


png('/home2/sbjs0428/R/atopyNK/plots/batch/skin/cellchat/atsk_atsktx_outandin_heatmap_11x16.png', width=11, height=16, units='in', res=1200)
p1 <- netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = "all", signaling = pathway.union, title = names(object.list)[1], width = 10, height = 33, color.use=unname(skin_colors[levels(object.list[[1]]@meta$labels)]))
p2 <- netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = "all", signaling = pathway.union, title = names(object.list)[2], width = 10, height = 33, color.use=unname(skin_colors[levels(object.list[[2]]@meta$labels)]))
draw(p1+p2, ht_gap=unit(0.5, 'cm'))
dev.off()



#draw bubble plot for all cells 

for (cell in levels(cellchat@meta$labels)) {
  index <- which(levels(cellchat@meta$labels) == cell)
  dir_path <- paste0('/home2/sbjs0428/R/atopyNK/plots/batch/skin/cellchat/bubble/', cell)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  # Handle increased signaling plot
  tryCatch({
    plot <- netVisual_bubble(cellchat, sources.use = index, comparison = c(1, 2),
                             remove.isolate = TRUE, angle.x = 45, max.dataset = 2,
                             title.name = 'Increased Signaling in ATSKTX')
    ggsave(paste0(dir_path, '/atsk_atsktx_bubble_', cell, '_increased_10x18.png'),
           plot, width = 10, height = 18, dpi = 1200, units = 'in')
  }, error = function(e) {
    message(paste("Error in increased signaling plot for cell:", cell))
    message(e)
    unlink(dir_path, recursive = TRUE)
  })

  # Handle decreased signaling plot
  tryCatch({
    plot <- netVisual_bubble(cellchat, sources.use = index, comparison = c(1, 2),
                             remove.isolate = TRUE, angle.x = 45, max.dataset = 1,
                             title.name = 'Decreased Signaling in ATSKTX')
    ggsave(paste0(dir_path, '/atsk_atsktx_bubble_', cell, '_decreased_10x18.png'),
           plot, width = 10, height = 18, dpi = 1200, units = 'in')
  }, error = function(e) {
    message(paste("Error in decreased signaling plot for cell:", cell))
    message(e)
    unlink(dir_path, recursive = TRUE)
  })
}




#### three samples 

objlist <- readRDS('/home2/sbjs0428/R/atopyNK/cellchat/batch/skin/data/objlist.RDS')
n <- readRDS('/home2/sbjs0428/R/atopyNK/cellchat/batch/skin/data/NSK_cellchat.RDS')
oblist <- list('NSK' = n, "ATSK" = objlist[[1]], 'ATSKTX' = objlist[[2]])
cellchat <- mergeCellChat(oblist, add.names= names(oblist))

png('/home2/sbjs0428/R/atopyNK/plots/batch/skin/cellchat/nsk_atsk_atsktx_pathway_barplot_10x18.png', width=10, height=18, units='in', res=1200)
rankNet(cellchat, mode='comparison', measure='weight', sources.use=NULL, targets.use=NULL, stacked=T, do.stat=T, color.use=c('#00BE67',"#377EB8", "#E41A1C"), comparison=c(1,2,3))
dev.off()

library(ggplot2)
library(cowplot)

cellchat <- readRDS('/home2/sbjs0428/R/atopyNK/cellchat/batch/skin/data/cellchat_merged_at_atx.RDS')

pathways = c('PSAP', 'FGF', 'NRG', 'NKG2D', 'GDF', 'PVR', 'IFN-II', 'CHEMERIN', 'LIFR', 'NT', '2-AG', 'EPGN', 'NEGR', 'CSPG4', 'TGFb', 'TWEAK', 'CSF', 'OCLN', 'IL17', 'VISFATIN', 'CD96', 'CD39', 'GAP', 'IL6', 'THBS', 'ANGPTL', 'WNT', 'VEGF', 'IL4', 'COLLAGEN', 'CCL', 'KLK', 'CXCL', 'PECAM2', 'MIF')
#pathways = c('PSAP', 'FGF', 'GDF', 'CHEMERIN', 'PDL2', 'PVR', 'PTPR', '12oxoLTB5', '2-AG', 'NEGR', 'Adenosine', 'CD39', 'NT', 'TGFb', 'IL6', 'IL17', 'WNT', 'PECAM2', 'EGF', 'CXCL', 'WNT', 'NRG', 'CCL', 'IFN-II')
sam = 'ATSK'

for (pathway in pathways) {
  path = '/home2/sbjs0428/R/atopyNK/plots/batch/skin/cellchat/vlnplot/atsk/'
  if (!dir.exists(path)) {
	dir.create(path, recursive=T)
  }
  # Set up the PNG device with the specified filename and dimensions
  png(paste0(path, pathway, '_', sam, '_enriched', '_vlnplot_12x12.png'), 
      width = 12, height = 12, units = 'in', res = 1200)
  
  # Generate the plot
  plot <- plotGeneExpression(cellchat, signaling = pathway, split.by = 'datasets', 
                             type = 'violin', color.use = c("#377EB8", "#E41A1C"), angle.x = 45)
  
  # Add title and adjust margins
  plot <- plot_grid(plot, ncol = 1, align = 'v') + 
          ggtitle(paste0(pathway, ' pathway genes ATSK vs ATSKTX')) + 
          theme(plot.title = element_text(size = 15, face = 'bold', hjust = 0.5), 
                plot.margin = margin(t = 0.2, r = 0, b = 0.13, l = 0, unit = "in"))
  
  # Print the plot
  print(plot)
  
  # Close the PNG device
  dev.off()
}

pathways <- c('Cholesterol', 'CDH1', 'HSPG', 'THY1', 'PD-L1', 'NPR1', 'CDH', 'KIT', 'SPP1', 'SEMA4', 'NOTCH', 'Netrin', 'DHT', 'EPHB', 'ADGRL', 'IGF', 'BMP', 'SLIT', 'CypA', 'UNC5', 'HH', 'VCAM', 'DHEA', 'MAG', 'CEACAM', 'BST2', 'CD34', 'NRXN', 'SELPLG', 'MHC-II', 'SELE')


#pathways = c('IGF', 'CypA', 'Prostaglandin', 'SLURP', 'FLT3', 'HH', 'DHEA', 'CD40', 'MAG', 'VCAM', 'CEACAM', 'BST2', 'CD34', 'TNF', 'NRXN', 'SELPLG', 'MHC-II', 'SELE', 'NOTCH', 'BMP')

sam = 'ATSKTX'
for (pathway in pathways) {
  # Set up the PNG device with the specified filename and dimensions
  path = '/home2/sbjs0428/R/atopyNK/plots/batch/skin/cellchat/vlnplot/atsktx/'
  if (!dir.exists(path)) {
	dir.create(path, recursive=T)
  }
  png(paste0(path, pathway, '_', sam, '_enriched', '_vlnplot_12x12.png'),
      width = 12, height = 12, units = 'in', res = 1200)

  # Generate the plot
  plot <- plotGeneExpression(cellchat, signaling = pathway, split.by = 'datasets',
                             type = 'violin', color.use = c("#377EB8", "#E41A1C"), angle.x = 45)

  # Add title and adjust margins
  plot <- plot_grid(plot, ncol = 1, align = 'v') +
          ggtitle(paste0(pathway, ' pathway genes ATSK vs ATSKTX')) +
          theme(plot.title = element_text(size = 15, face = 'bold', hjust = 0.5),
                plot.margin = margin(t = 0.2, r = 0, b = 0.13, l = 0, unit = "in"))

  # Print the plot
  print(plot)

  # Close the PNG device
  dev.off()
}


#####################


#cellchat <- identifyOverExpressedGenes(cellchat, group.by='datasets', pos.dataset=pos.dataset, features.name=features.name, only.pos=F, group.DE.combined=F)
#net <- netMappingDEG(cellchat, features.name=features.name, variable.all=T)
