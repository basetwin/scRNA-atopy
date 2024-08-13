source('~/necessary_packages.R')
obj <- readRDS('/home2/sbjs0428/R/atopyNK/workspace/nobatch/lymph_obj.RDS')

atln_deg <- FindMarkers(obj, group.by = 'Type', ident.1 = 'ATLN', ident.2 = 'ATLNTX')
atlntx_deg <- FindMarkers(obj, group.by = 'Type', ident.1 = 'ATLNTX', ident.2 = 'ATLN')

atln_deg_names <- rownames(atln_deg %>%  filter(p_val_adj < 0.05) %>% filter(avg_log2FC > 1))

atlntx_deg_names <- rownames(atlntx_deg %>%  filter(p_val_adj < 0.05) %>% filter(avg_log2FC > 1))

selected_genes <- c(atlntx_deg_names, atln_deg_names)
groups <- list(ATLNTX = 1:29, ATLN = 30:101)

# Extract expression data for the selected genes from the Seurat object
expr_data <- obj[['RNA']]$data[selected_genes, ]

# Convert to a data frame and transpose
expr_data_df <- as.data.frame(expr_data)
expr_data_t <- t(expr_data_df)

# Compute Spearman correlation matrix
CorrelationMatrix <- cor(expr_data_t, method = "spearman")

# Specify the file path and name

output_file_path <- '/home2/sbjs0428/R/atopyNK/plots/nobatch/lymph/lymph_type_deg_network_legend.png'

# Open a PNG device
png(filename = output_file_path, width = 15, height = 15, units = 'in', res = 1200)

# Plot
qgraph(CorrelationMatrix, graph = "cor", layout = "spring", vsize=4,
       directed = FALSE, posCol = "brown4", negCol = "blue4",
       edge.labels = F,
       labels=colnames(CorrelationMatrix),
       label.prop=1,
       minimum="sig", sampleSize=10559,  # Adjust sampleSize according to your data
       groups = groups,
#       color = c("#FF3366", "#99FF00"),
       color = c("#E41A1C","#377EB8"), 
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
#atln_deg <- atln_deg %>%  filter(p_val_adj < 0.05) %>% filter(avg_log2FC > 1)
#
#atlntx_deg <- atlntx_deg %>%  filter(p_val_adj < 0.05) %>% filter(avg_log2FC > 1)
#
#atln_deg$type <- 2
#atlntx_deg$type <- 1
#atln_deg$gene <- rownames(atln_deg)
#atlntx_deg$gene <- rownames(atlntx_deg)
#
#deg <- rbind(atlntx_deg, atln_deg)
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
##'#FF3366' -> red (atlntx)
##'#99FF00' -> green (atln)
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
##/home2/sbjs0428/R/atopyNK/plots/nobatch/lymph
#output_file_path <- '/home2/sbjs0428/R/atopyNK/plots/nobatch/lymph/lymph_type_deg_protein_network_10x10.png'
#
## Open a PNG device
#png(filename = output_file_path, width = 10, height = 10, units = 'in', res = 1000)
#
#
##plotting
#string_db$plot_network(hits, payload_id = payload_id)
#
#dev.off()

