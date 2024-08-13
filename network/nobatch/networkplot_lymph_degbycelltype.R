source('~/necessary_packages.R')
obj <- readRDS('/home2/sbjs0428/R/atopyNK/workspace/nobatch/lymph_obj.RDS')

cells <- list()
cell_type <- unique(obj$celltype)

for (type in cell_type) {
	cells[[type]] <- subset(obj, celltype == type)
}

cells_deg <- list()

for (cell_name in names(cells)) {
	filtered_deg <- FindMarkers(cells[[cell_name]], group.by='Type', ident.1='ATLNTX', ident.2='ATLN', only.pos=T) %>%
			filter(p_val_adj < 0.05, avg_log2FC > 1) %>% arrange(desc(avg_log2FC)) %>% slice_head(n=50)
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


color <- lymph_colors

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

new_colors[new_colors == "#E69F00"] <- "#cc7b02"

# Extract expression data for the selected genes from the Seurat object
expr_data <- obj[['RNA']]$counts[selected_genes, ]

# Convert to a data frame and transpose
expr_data_df <- as.data.frame(expr_data)
expr_data_t <- t(expr_data_df)

# Compute Spearman correlation matrix
CorrelationMatrix <- cor(expr_data_t, method = "spearman")

# Specify the file path and name
#output_file_path <- '/home2/sbjs0428/R/atopyNK/plots/lymph/lymph_type_degpercelltype_network.png'
output_file_path <- 'temp1.png'


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



#library(STRINGdb)
#string_db <- STRINGdb$new(version="11.5", species=10090,score_threshold=200, input_directory="")
#options(download.file.method="libcurl")
#
#deg <- cells_deg
#
#celltypes <- unique(deg$celltype)
#
#for (type in celltypes) {
#	deg[deg$celltype == type, 'celltype'] <- which(celltypes == type)
#}
#
#deg$celltype <- as.numeric(as.character(deg$celltype))
#deg$gene <- rownames(deg)
#
#example1_mapped <- string_db$map(deg, "gene", removeUnmappedRows = TRUE)
#hits <- example1_mapped$STRING_id
#
#example1_mapped_logFC <- string_db$add_diff_exp_color(example1_mapped, logFcColStr="celltype")
#unique(example1_mapped_logFC$color)
#
#for (type in celltypes) {
#	example1_mapped_logFC$color2[example1_mapped_logFC$celltype == which(celltypes == type)] <-  color[[type]]
#}
#payload_id <- string_db$post_payload(example1_mapped_logFC$STRING_id, 
#                                     colors = example1_mapped_logFC$color2)
#
##output_file_path <- '/home2/sbjs0428/R/atopyNK/plots/lymph/lymph_type_degpercelltype_protein_network_10x10.png'
#output_file_patj <- 'temp_2.png'
#
#png(filename = output_file_path, width = 10, height = 10, units = 'in', res = 1638)
#
##plotting
#string_db$plot_network(hits, payload_id = payload_id)
#dev.off()
