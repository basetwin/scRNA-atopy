source('~/necessary_packages.R')
obj <- readRDS('/home2/sbjs0428/R/atopyNK/workspace/nobatch/skin_obj.RDS')

###cellchat new color
#E41A1C <- atlntx (red)  "#FFFFFFFF" 
#377EB8 <- atln ( blue)  "#FF0000FF"


atsk_deg <- FindMarkers(obj, group.by = 'Type', ident.1 = 'ATSK', ident.2 = 'ATSKTX')
atsktx_deg <- FindMarkers(obj, group.by = 'Type', ident.1 = 'ATSKTX', ident.2 = 'ATSK')

atsk_deg_names <- rownames(atsk_deg %>%  filter(p_val_adj < 0.001) %>% filter(avg_log2FC > 1.5))

atsktx_deg_names <- rownames(atsktx_deg %>%  filter(p_val_adj < 0.001) %>% filter(avg_log2FC > 1.5))

selected_genes <- c(atsktx_deg_names, atsk_deg_names)
groups <- list(ATSKTX = 1:52, ATSK = 53:120)

# Extract expression data for the selected genes from the Seurat object
expr_data <- obj[['RNA']]$data[selected_genes, ]

# Convert to a data frame and transpose
expr_data_df <- as.data.frame(expr_data)
expr_data_t <- t(expr_data_df)

# Compute Spearman correlation matrix
CorrelationMatrix <- cor(expr_data_t, method = "spearman")

# Specify the file path and name
#output_file_path <- '/home2/sbjs0428/R/atopyNK/plots/nobatch/skin/skin_type_deg_network.png'
output_file_path <- 'skin_legend.png'

# Open a PNG device
png(filename = output_file_path, width = 15, height = 15, units = 'in', res = 1000)

# Plot
qgraph(CorrelationMatrix, graph = "cor", layout = "spring", vsize=4,
       directed = FALSE, posCol = "brown4", negCol = "blue4",
       edge.labels = F,
       labels=colnames(CorrelationMatrix),
       label.prop=1,
       minimum="sig", sampleSize=10559,  # Adjust sampleSize according to your data
       groups = groups,
       color = c("#E41A1C", "#377EB8"),
       border.color = "white", border.width = 1,
       edge.width = 0.8,
       curve = 0.2, curveAll = T,
       label.color = "black",
       legend = T,
       label.cex = 1.5,
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
#atsk_deg <- atsk_deg %>%  filter(p_val_adj < 0.001) %>% filter(avg_log2FC > 1.5)
#
#atsktx_deg <- atsktx_deg %>%  filter(p_val_adj < 0.001) %>% filter(avg_log2FC > 1.5)
#
#atsk_deg$type <- 2
#atsktx_deg$type <- 1
#
#deg <- rbind(atsktx_deg, atsk_deg)
#deg$type <- as.numeric(as.character(deg$type))
#deg$gene <- rownames(deg)
#
#
#example1_mapped <- string_db$map(deg, "gene", removeUnmappedRows = TRUE)
#hits <- example1_mapped$STRING_id
#
#example1_mapped_logFC <- string_db$add_diff_exp_color(example1_mapped, logFcColStr="type")
#unique(example1_mapped_logFC$color)
##여기서 color에 따라 지정 컬러 바꿔주자
#
##'#FF3366' -> red (atsktx)
##'#99FF00' -> green (atsk)
#
####cellchat new color
##E41A1C <- atlntx (red)  "#FFFFFFFF" 
##377EB8 <- atln ( blue)  "#FF0000FF"
#
#
#
#
#example1_mapped_logFC$color2[example1_mapped_logFC$color == "#FFFFFFFF"] <- '#E41A1C'
#example1_mapped_logFC$color2[example1_mapped_logFC$color == "#FF0000FF"] <- '#377EB8'
#payload_id <- string_db$post_payload(example1_mapped_logFC$STRING_id, 
#                                     colors = example1_mapped_logFC$color2)
#
#output_file_path <- '/home2/sbjs0428/R/atopyNK/plots/nobatch/skin/skin_type_deg_protein_network_10x10.png'
#png(filename = output_file_path, width = 10, height = 10, units = 'in', res = 1200)
#
##plotting
#string_db$plot_network(hits, payload_id = payload_id)
#dev.off()
