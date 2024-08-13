source('~/necessary_packages.R')
load('skin_240318.RData')

cells <- list()
cell_type <- unique(obj$celltype)

cells[['Keratinocyte']] <- subset(obj, celltype %in% c(as.character(cell_type[grepl('keratinocyte', cell_type)]), 'Bulge cell'))
cells[['Endothelial cell']] <- subset(obj, celltype %in% cell_type[grepl('Endothelial', cell_type)])
cells[['Th cell']] <- subset(obj, celltype %in% cell_type[grepl('Th', cell_type)])
cells[['Adipocyte']] <- subset(obj, celltype == 'Adipocyte')
cells[['Langerhans cell']] <- subset(obj, celltype == 'Langerhans cell')
cells[['Melanocyte']] <- subset(obj, celltype == 'Melanocyte')
cells[['Fibroblast']] <- subset(obj, celltype == 'Fibroblast')
cells[['Schwan cell']] <- subset(obj, celltype == 'Schwan cell')
cells[['NK cell']] <- subset(obj, celltype == 'NK cell')
cells[['Sebaceous gland']] <- subset(obj, celltype == 'Sebaceous gland')

cells_deg <- list()

for (cell_name in names(cells)) {
	filtered_deg <- FindMarkers(cells[[cell_name]], group.by='Type', ident.1='ATSKTX', ident.2='ATSK', only.pos=T) %>%
			filter(p_val_adj < 0.05) %>% filter(avg_log2FC > 0.6) %>% arrange(desc(pct.1-pct.2)) %>% slice_head(n=50)
	if (nrow(filtered_deg) > 0) {
		cells_deg[[cell_name]] <- filtered_deg
		cells_deg[[cell_name]]$celltype <- cell_name 
		cells_deg[[cell_name]]$gene <- rownames(cells_deg[[cell_name]])
	} 
}

cells_deg <- do.call(rbind, cells_deg)

cells_deg <- cells_deg %>% group_by(gene) %>% arrange(desc(avg_log2FC)) %>% slice_head(n=1) %>% ungroup()
cells_deg <- as.data.frame(cells_deg)
rownames(cells_deg) <- cells_deg$gene
cells_deg$gene <- NULL

 
#colors by cell type
#
#                  Bulge cell              Sebaceous gland 
#                   "#E69F00"                    "#56B4E9" 
#       Krt6a- B keratinocyte Krt6a- Krt79+ B keratinocyte 
#                   "#009E73"                    "#F0E442" 
#       Krt6a- S keratinocyte       Krt79+ BS keratinocyte 
#                   "#0072B2"                    "#D55E00" 
#              B keratinocyte              BS keratinocyte 
#                   "#CC79A7"                    "#666666" 
#Proliferating B keratinocyte               G keratinocyte 
#                   "#AD7700"                    "#1C91D4" 
#                   Adipocyte     Ccl21a+ Endothelial cell 
#                   "#007756"                    "#D5C711" 
#   Cd300lg+ Endothelial cell                   Fibroblast 
#                   "#005685"                    "#A04700" 
#             Langerhans cell                   Melanocyte 
#                   "#B14380"                    "#4D4D4D" 
#                 Schwan cell                          Th2 
#                   "#FFBE2D"                    "#80C7EF" 
#                        Th17           Proliferating Th17 
#                   "#00F6B3"                    "#F4EB71" 
#                        Th22                      NK cell 
#                   "#06A5FF"                    "#FF8320" 


color <- dittoHeatmap(obj, c('Flg', 'Lor', 'Ivl'), annot.by = 'celltype', data.out=T)$annotation_colors$celltype
color['Keratinocyte'] <- "#CC79A7"
color['Endothelial cell'] <- "#005685"
color['Th cell'] <- "#F4EB71"

selected_genes <- rownames(cells_deg)	
groups <- list()
for (type in unique(cells_deg$celltype)) {
	index <- which(cells_deg$celltype == type)
	groups[[type]] <- index
}
new_colors <- vector()

for (type in names(groups)) {
	new_colors <- c(new_colors, color[[type]]) 
}

# Extract expression data for the selected genes from the Seurat object
expr_data <- obj[['RNA']]$counts[selected_genes, ]

# Convert to a data frame and transpose
expr_data_df <- as.data.frame(expr_data)
expr_data_t <- t(expr_data_df)

# Compute Spearman correlation matrix
CorrelationMatrix <- cor(expr_data_t, method = "spearman")

# Specify the file path and name
output_file_path <- '/home2/sbjs0428/R/atopyNK/plots/skin/skin_type_degpercelltype_network.png'

# Open a PNG device
png(filename = output_file_path, width = 20, height = 15, units = 'in', res = 1638)

# Plot
qgraph(CorrelationMatrix, graph = "cor", layout = "spring", vsize=4,
       directed = FALSE, posCol = "brown4", negCol = "blue4",
       edge.labels = F,
       labels=colnames(CorrelationMatrix),
       label.prop=1,
       minimum="sig", sampleSize=10559,  # Adjust sampleSize according to your data
       groups = groups,
       color = new_colors,
       border.color = "white", border.width = 1,
       edge.width = 0.8,
       curve = 0.2, curveAll = T,
       label.color = "black",
       legend = T,
       label.cex = 2,
       label.scale.equal = T,
       repulsion=2,
       node.width=0.5, node.height=0.5
       )

# Close the device
dev.off()



library(STRINGdb)
string_db <- STRINGdb$new(version="11.5", species=10090,score_threshold=200, input_directory="")
options(download.file.method="libcurl")

deg <- cells_deg

# celltype order
#[1] "NK cell"          "Adipocyte"        "Keratinocyte"     "Fibroblast"      
#[5] "Th cell"          "Melanocyte"       "Sebaceous gland"  "Endothelial cell"

celltypes <- unique(deg$celltype)

for (type in celltypes) {
	deg[deg$celltype == type, 'celltype'] <- which(celltypes == type)
}

deg$celltype <- as.numeric(as.character(deg$celltype))
deg$gene <- rownames(deg)

example1_mapped <- string_db$map(deg, "gene", removeUnmappedRows = TRUE)
hits <- example1_mapped$STRING_id

example1_mapped_logFC <- string_db$add_diff_exp_color(example1_mapped, logFcColStr="celltype")
unique(example1_mapped_logFC$color)

for (type in celltypes) {
	example1_mapped_logFC$color2[example1_mapped_logFC$celltype == which(celltypes == type)] <-  color[[type]]
}
payload_id <- string_db$post_payload(example1_mapped_logFC$STRING_id, 
                                     colors = example1_mapped_logFC$color2)

output_file_path <- '/home2/sbjs0428/R/atopyNK/plots/skin/skin_type_degpercelltype_protein_network_10x10.png'
png(filename = output_file_path, width = 10, height = 10, units = 'in', res = 1638)

#plotting
string_db$plot_network(hits, payload_id = payload_id)
dev.off()
