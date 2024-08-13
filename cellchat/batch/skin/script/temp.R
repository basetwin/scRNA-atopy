source('~/necessary_packages.R')

obj <- readRDS("/home2/sbjs0428/R/atopyNK/workspace/nobatch/skin_obj_240708.RDS")


#DimPlot drawing

skin_colors <- c(
  "Bulge cell" = "#E69F00",
  "Sebaceous gland" = "#56B4E9",
  "Krt6a- B keratinocyte" = "#009E73",
  "Krt6a- Krt79+ B keratinocyte" = "#F0E442",
  "Krt6a- BS keratinocyte" = "#0072B2",
  "Krt6a- S keratinocyte" = "#D55E00",
  "Hyperproliferative B keratinocyte" = "#CC79A7",
  "Krt79+ BS keratinocyte" = "#666666",
  "BS keratinocyte" = "#AD7700",
  "G keratinocyte" = "#1C91D4",
  "Adipocyte" = "#007756",
  "Ccl21a+ Endothelial cell" = "#D5C711",
  "Cd300lg+ Endothelial cell" = "#005685",
  "Fibroblast" = "#A04700",
  "Langerhans cell" = "#B14380",
  "Melanocyte" = "#4D4D4D",
  "Schwan cell" = "#FFBE2D",
  "Treg" = "#80C7EF",
  "Th17" = "#00F6B3",
  "Proliferating Th17" = "#F4EB71",
  "Th22" = "#06A5FF",
  "NK cell" = "#FF8320"
)


#object.list <- readRDS('/home2/sbjs0428/R/atopyNK/cellchat/nobatch/skin/data/objlist.RDS')
#cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#
#reticulate::use_python('~/miniconda3/envs/monocle3/bin/python', required=T)
#
##cellchat <- netAnalysis_computeCentrality(cellchat, slot.name='netP')
#
#cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
#
#cellchat <- netEmbedding(cellchat, type = "functional")
#cellchat <- netClustering(cellchat, type = "functional")
#cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
#cellchat <- netEmbedding(cellchat, type = "structural")
#cellchat <- netClustering(cellchat, type = "structural")
#
#
#
#cellchat@idents$joint <- factor(cellchat@idents$joint, levels= levels(obj$celltype))
#
#saveRDS(cellchat, '/home2/sbjs0428/R/atopyNK/cellchat/nobatch/skin/data/cellchat_merged_at_atx.RDS')
#
#pathway.union <- union(object.list[[1]]@netP$pathways, object.list[[2]]@netP$pathways)
#
#
#
#library(ComplexHeatmap)
#png('/home2/sbjs0428/R/atopyNK/plots/nobatch/skin/cellchat/atsk_atsktx_pathway_barplot_8x18.png', width=8, height=18, units='in', res=1200)
#rankNet(cellchat, mode='comparison', measure='weight', sources.use=NULL, targets.use=NULL, stacked=T, do.stat=T, color.use=c("#377EB8", "#E41A1C"))
#dev.off()
#
#
#source('/home2/sbjs0428/R/atopyNK/cellchat/nobatch/skin/script/netVisual_heatmap.R')
#png('/home2/sbjs0428/R/atopyNK/plots/nobatch/skin/cellchat/atsk_atsktx_differential_weight_heatmap_12x10.png', width=12, height=10, units='in', res=1200)
#netVisual_heatmap(cellchat, measure='weight', title.name='Differential interactions-weight ATSK vs ATSKTX', color.use=skin_colors)
#dev.off()
#
#png('/home2/sbjs0428/R/atopyNK/plots/nobatch/skin/cellchat/atsk_atsktx_differential_count_heatmap_12x10.png', width=12, height=10, units='in', res=1200)
#netVisual_heatmap(cellchat, measure='count', title.name='Differential interactions-count ATSK vs ATSKTX', color.use=skin_colors)
#dev.off()
#
#png('/home2/sbjs0428/R/atopyNK/plots/nobatch/skin/cellchat/atsk_atsktx_outgoing_heatmap_11x16.png', width=11, height=16, units='in', res=1200)
#
#p1 <- netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[1], width = 10, height = 33, color.use=unname(skin_colors[levels(object.list[[1]]@meta$labels)]))
#
#p2 <- netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[2], width = 10, height = 33, color.use=unname(skin_colors[levels(object.list[[2]]@meta$labels)]))
#
#draw(p1+p2, ht_gap=unit(0.5, 'cm'))
#dev.off()

cellchat <- readRDS('/home2/sbjs0428/R/atopyNK/cellchat/nobatch/skin/data/cellchat_merged_at_atx.RDS')
object.list <- readRDS('/home2/sbjs0428/R/atopyNK/cellchat/nobatch/skin/data/objlist.RDS')
obj <- readRDS('/home2/sbjs0428/R/atopyNK/workspace/nobatch/skin_obj_240708.RDS')
pathway.union <- union(object.list[[1]]@netP$pathways, object.list[[2]]@netP$pathways)
library(ComplexHeatmap)
source('/home2/sbjs0428/R/atopyNK/cellchat/nobatch/skin/script/netVisual_heatmap.R')




png('/home2/sbjs0428/R/atopyNK/plots/nobatch/skin/cellchat/atsk_atsktx_incoming_heatmap_11x16.png', width=11, height=16, units='in', res=1200)
p1 <- netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[1], width = 10, height = 33, color.use=unname(skin_colors[levels(object.list[[1]]@meta$labels)]))
p2 <- netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[2], width = 10, height = 33, color.use=unname(skin_colors[levels(object.list[[2]]@meta$labels)]))
draw(p1+p2, ht_gap=unit(0.5, 'cm'))
dev.off()


png('/home2/sbjs0428/R/atopyNK/plots/nobatch/skin/cellchat/atsk_atsktx_outandin_heatmap_11x16.png', width=11, height=16, units='in', res=1200)
p1 <- netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = "all", signaling = pathway.union, title = names(object.list)[1], width = 10, height = 33, color.use=unname(skin_colors[levels(object.list[[1]]@meta$labels)]))
p2 <- netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = "all", signaling = pathway.union, title = names(object.list)[2], width = 10, height = 33, color.use=unname(skin_colors[levels(object.list[[2]]@meta$labels)]))
draw(p1+p2, ht_gap=unit(0.5, 'cm'))
dev.off()

#draw bubble plot for all cells 

for (cell in levels(cellchat@meta$labels)) {
  index <- which(levels(cellchat@meta$labels) == cell)
  dir_path <- paste0('/home2/sbjs0428/R/atopyNK/plots/nobatch/skin/cellchat/bubble/', cell)
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

objlist <- readRDS('/home2/sbjs0428/R/atopyNK/cellchat/nobatch/skin/data/objlist.RDS')
n <- readRDS('/home2/sbjs0428/R/atopyNK/cellchat/nobatch/skin/data/NSK_cellchat.RDS')
oblist <- list('NSK' = n, "ATSK" = objlist[[1]], 'ATSKTX' = objlist[[2]])
cellchat <- mergeCellChat(oblist, add.names= names(oblist))

png('/home2/sbjs0428/R/atopyNK/plots/nobatch/skin/cellchat/nsk_atsk_atsktx_pathway_barplot_10x18.png', width=10, height=18, units='in', res=1200)
rankNet(cellchat, mode='comparison', measure='weight', sources.use=NULL, targets.use=NULL, stacked=T, do.stat=T, color.use=c('#00BE67',"#377EB8", "#E41A1C"), comparison=c(1,2,3))
dev.off()
##################


library(ggplot2)
library(cowplot)

cellchat <- readRDS('/home2/sbjs0428/R/atopyNK/cellchat/nobatch/skin/data/cellchat_merged_at_atx.RDS')

pathways = c('PSAP', 'FGF', 'GDF', 'CHEMERIN', 'PDL2', 'PVR', 'PTPR', '12oxoLTB5', '2-AG', 'NEGR', 'Adenosine', 'CD39', 'NT', 'TGFb', 'IL6', 'IL17', 'WNT', 'PECAM2', 'EGF', 'CXCL', 'WNT', 'NRG')
sam = 'ATSK'

for (pathway in pathways) {
  # Set up the PNG device with the specified filename and dimensions
  png(paste0('/home2/sbjs0428/R/atopyNK/plots/nobatch/skin/cellchat/vlnplot/atsk/', pathway, '_', sam, '_enriched', '_vlnplot_12x12.png'),
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

pathways = c('IGF', 'CypA', 'Prostaglandin', 'SLURP', 'FLT3', 'HH', 'DHEA', 'CD40', 'MAG', 'VCAM', 'CEACAM', 'BST2', 'CD34', 'TNF', 'NRXN', 'SELPLG', 'MHC-II', 'SELE', 'NOTCH', 'BMP')
sam = 'ATSKTX'
for (pathway in pathways) {
  # Set up the PNG device with the specified filename and dimensions
  png(paste0('/home2/sbjs0428/R/atopyNK/plots/nobatch/skin/cellchat/vlnplot/atsktx/', pathway, '_', sam, '_enriched', '_vlnplot_12x12.png'),
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




