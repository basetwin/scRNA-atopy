library(Seurat)
library(Matrix)
library(tidyverse)
library(scDblFinder)
#library(future)
library(pheatmap)
library(openxlsx)
#library(BiocParallel)
library(ggplot2)
library(dplyr)
library(dittoSeq)
library(clusterProfiler)
library(enrichplot)

dirs <-  list.dirs(path = '/home2/sbjs0428/R/atopyNK/copy/lymph', recursive = F, full.names = F)

for(x in dirs){
  name <- x 
  #name <- gsub('_filtered_feature_bc_matrix','', x) 
  # x -> characters in dirs / _filtered_feature_bc_matrix -> expression pattern that the 
  #'gsub' function is looking for in the string 'x'/ ' ' -> is the replacement string.  
  # if a pattern '_filtered_feature_bc_matrix' is found in 'x', it will be replaced with a space ' '
  
  
  cts <- ReadMtx(mtx = paste0('/home2/sbjs0428/R/atopyNK/copy/lymph/',x,'/matrix.mtx.gz'),
                 features = paste0('/home2/sbjs0428/R/atopyNK/copy/lymph/',x,'/features.tsv.gz'),
                 cells = paste0('/home2/sbjs0428/R/atopyNK/copy/lymph/',x,'/barcodes.tsv.gz'))
  # ReadMtx is a function used to read data from the matrix, features, and barcodes files
  # mtx, features, and cells arguments specify the paths to these files.
  # The data read from these files is stored in the variable cts.
  # assigns this newly created Seurat object to a variable with the name stored in the name variable.
  # Essentially, it gives a meaningful name to each Seurat object based on the directory name 
  obj <- CreateSeuratObject(counts = cts)
  sce <- scDblFinder(obj[['RNA']]$counts)
  obj$doublet <- sce$scDblFinder.class
  obj <- subset(obj, subset = doublet == 'singlet')
  assign(name, obj)
  
}

merged_seurat <- merge(normallymphnode, y = c(atopylymphnode, atopylymphnodetx),
                       add.cell.ids = c('NLN','ATLN', 'ATLNTX'),
                       project = 'AL')
# create a sample column
merged_seurat$sample <- rownames(merged_seurat@meta.data)

merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = 'sample', into = c('Type', 'Barcode'), 
                                    sep = '_')
  


# calculate mitochondrial percentage
merged_seurat$mitoPercent <- PercentageFeatureSet(merged_seurat, pattern='^mt-')
merged_seurat_filtered <- subset(merged_seurat, subset = nCount_RNA > 500 &
                                   nFeature_RNA > 300 &
                                   mitoPercent < 10)

rm(atopylymphnode, atopylymphnodetx, normallymphnode, name, x)

#DimPlot(merged_seurat_filtered, reduction = 'umap')
rm(dirs, sce,cts)

# Calculate the number of cells in which each gene is expressed
gene_counts <- rowSums(merged_seurat[['RNA']]$counts > 0)

# Get genes that are expressed in at least 3 cells
genes_above_threshold <- names(gene_counts[gene_counts >= 3])


# Combine the genes above threshold with the genes in gene_list
final_gene_list <- genes_above_threshold

# Subset the Seurat object to retain only the genes in final_gene_list
merged_seurat <- subset(merged_seurat, features = final_gene_list)
rm(final_gene_list, genes_above_threshold, gene_counts)


# Extract all gene names
allgenes <- rownames(merged_seurat[['RNA']]$counts)

# Use grepl to create a logical vector for each condition
starts_with_Gm_or_Mt <- grepl("^(Gm|Mt|mt)", allgenes)
ends_with_Rik <- grepl("Rik$", allgenes)

# Combine conditions: genes that start with 'Gm' or 'Mt' OR end with 'Rik'
genes_to_exclude <- starts_with_Gm_or_Mt | ends_with_Rik

# Subset the allgenes vector to exclude these genes
filtered_genes <- allgenes[!genes_to_exclude]

allgenes_filtered <- filtered_genes
rm(starts_with_Gm_or_Mt, ends_with_Rik, genes_to_exclude, filtered_genes, allgenes)
merged_seurat <- subset(merged_seurat, features = allgenes_filtered)
rm(allgenes_filtered)


merged_seurat <- NormalizeData(merged_seurat)
merged_seurat <- FindVariableFeatures(object = merged_seurat)
merged_seurat <- ScaleData(object = merged_seurat, features = rownames(merged_seurat))
merged_seurat <- RunPCA(object = merged_seurat)
merged_seurat <- IntegrateLayers(merged_seurat, method =RPCAIntegration, orig.reduction='pca', new.reduction='integrated_rpca', verbose=F)
merged_seurat <- FindNeighbors(object = merged_seurat, dims = 1:20, reduction='integrated_rpca')
merged_seurat <- FindClusters(object = merged_seurat, resolution = 0.8, verbose = FALSE, cluster.name = 'rpca_clusters')
merged_seurat <- RunUMAP(object = merged_seurat, dims = 1:20, reduction = 'integrated_rpca', reduction.name = 'umap_rpca')
merged_seurat <- JoinLayers(merged_seurat)



#############################################################################################
#############################################################################################


markers <- FindAllMarkers(merged_seurat_genefiltered, only.pos = TRUE, verbose = FALSE)


#filtering genes starting with 'Mt', 'Gm' and ending with 'Rik'
markers <- markers%>%
  filter(!grepl("^Gm", rownames(markers)) & 
           !grepl("Rik$", rownames(markers)) & 
           !grepl("^Mt", rownames(markers)))




# Assuming 'markers' is your marker table from FindAllMarkers
marklist <- list()

for (i in 0:22) {
  templist <- markers[markers$cluster == i, ]
  templist <- templist[order(templist$avg_log2FC, decreasing = TRUE), ]
  marklist[[i+1]] <- templist
}

rm(i, templist)


library(dplyr)

### The output genes from this process are not all in the data table  
# Assuming 'markers' is your marker table from FindAllMarkers

# marklist_top10 <- list()
# 
# for (i in 0:22) {
#   # Subset for each cluster
#   templist <- markers[markers$cluster == i, ]
#   
#   # Apply cutoff (e.g., avg_log2FC > 1)
#   templist <- templist %>% filter(avg_log2FC > 1)
#   
#   # Sort and select top 10 by avg_log2FC
#   templist_top10 <- templist %>% top_n(n = 10, wt = avg_log2FC)
#   
#   # Assign to list
#   marklist_top10[[i+1]] <- templist_top10
# }
# 
# top_10genes<- unique(unlist(lapply(marklist_top10, function(x) rownames(x))))
# rm(i, templist, templist_top10)



top30 <- markers %>%
  group_by(cluster) %>%
  filter(p_val_adj <= 0.05) %>%
  arrange(desc(abs(pct.1 - pct.2))) %>%
  slice_head(n = 30) %>%
  ungroup()


# Removing genes with names starting with 'Gm'
top30 <- top30 %>%
  filter(!grepl("^Gm", gene) & !grepl('Rik$', gene) & !grepl('^Mt', gene))


filtered_markers_between_clusters <- markers %>%
  group_by(cluster) %>%
  filter(p_val_adj <= 0.05) %>%
  filter(abs(pct.1 - pct.2) >= 0.8) %>%
  arrange(desc(abs(pct.1 - pct.2))) %>%
  ungroup()


#Find markers between conditions
markers_between_sample <- FindMarkers(merged_seurat_genefiltered, group.by = 'Type', ident.1 = 'ATLNTX', ident.2 = 'ATLN')
markers_between_sample <- markers_between_sample %>% arrange(desc(abs(pct.1 - pct.2)))

# Removing genes with names starting with 'Gm', 'Mt' and ending with 'Rik'
markers_between_sample <- markers_between_sample %>%
  filter(!grepl("^Gm", rownames(markers_between_sample)) & 
           !grepl("Rik$", rownames(markers_between_sample)) & 
           !grepl("^Mt", rownames(markers_between_sample)))


markers_between_sample <- markers_between_sample %>% filter(p_val_adj <= 0.05)


filtered_markers_between_types <- markers_between_sample %>% filter(abs(pct.1 - pct.2) >= 0.05)


#for filter list by expression level
filtered_markers_between_types_negative <- filtered_markers_between_types %>% filter(avg_log2FC <0) %>% filter(abs(pct.1 - pct.2) >= 0.1)

# intersection between top50 markers among clusters and markers of conditions
# cellmarkers <- rownames(markers)
# typemarkers <- rownames(markers_between_sample)
# cell_type_markers <- intersect(cellmarkers, typemarkers)
# rm(cellmarkers, typemarkers)


# get gene names to label in heatmap
set.seed(123)  # for reproducibility

# Assuming 'top30' is a dataframe with columns 'gene' and 'cluster'
sampled_genes <- top30 %>%
  group_by(cluster) %>%
  sample_n(size = 2) %>%
  ungroup()

# Extracting the gene names for the heatmap
genes_for_heatmap <- sampled_genes$gene

#drawing heatmap with labeling random gengenees
dittoHeatmap(
  object = merged_seurat_genefiltered,
  genes = celltype_markers_combined,  
  annot.by = c('celltype_cl', 'Type'),
  order.by = 'celltype_cl',
  show_rownames = F,
  scaled.to.max = TRUE,
)



#drawing heatmap using markers between conditions
dittoHeatmap(merged_seurat_genefiltered, rownames(markers_between_sample), annot.by = c('Type', 'celltype_cl'), order.by = 'celltype_cl', scaled.to.max = TRUE, show_rownames = FALSE)

#drawing heatmap using markers among seurat clusters 
dittoHeatmap(merged_seurat_genefiltered, top30$gene, annot.by = c('celltype_cl', 'Type'), order.by = 'celltype_cl', show_rownames = FALSE, scaled.to.max = TRUE)



#gene ontology

organism = "org.Mm.eg.db"
library(organism, character.only = TRUE)

# we want the log2 fold change 
original_gene_list <- filtered_markers_between_types_negative$avg_log2FC

# name the vector
names(original_gene_list) <- rownames(filteqred_markers_between_types_negative)

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# # List of genes to contain
# genes_to_exclude <- c("Ly6a", "Erlec1", "Jchain", "Ly6c2", "Ly6c1", "Lghg2b", "Lghg2c", "Cxcr3")
# 
# # Filter out the genes to exclude
# gene_list <- gene_list[names(original_gene_list) %in% genes_to_exclude]

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

rm(original_gene_list, genes_to_exclude)

gse <- gseGO(geneList=gene_list, 
             ont="ALL", 
             keyType="SYMBOL", 
             nPerm = 10000, 
             minGSSize=3, 
             maxGSSize=800, 
             pvalueCutoff=0.05, 
             verbose=TRUE, 
             OrgDb=organism, 
             pAdjustMethod="none")



require(DOSE)
dotplot(gse, split=".sign") + facet_grid(.~.sign)


data <- as.data.frame(merged_seurat_genefiltered@assays$RNA@data)
filtered_data <- data[top_50genes, ]

cellnames <- paste0(colnames(filtered_data), "_",merged_seurat_filtered@meta.data$celltype_cl)

colnames(filtered_data) <- cellnames
rm(cellnames)

sampletypes <- sapply(strsplit(colnames(filtered_data), "_"), `[`, 1)
celltypes <- sapply(strsplit(colnames(filtered_data), "_"), tail, n=1)
datacsv <- rbind(celltypes, filtered_data)
datacsv <- rbind(sampletypes, filtered_data)
rownames(filtered_data)[1:2] <- c("sampletype", "celltype")
rm(sampletypes, celltypes)


write.csv(filtered_data, "lymph_data.csv")

