library(Seurat)
library(Matrix)
library(tidyverse)
library(scDblFinder)
library(future)
library(pheatmap)
library(openxlsx)
library(BiocParallel)
library(ggpubr)
library(dittoSeq)
library(dplyr)
library(clusterProfiler)
library(enrichplot)
library(SingleR)



dirs <-  list.dirs(path = '/home2/sbjs0428/R/atopyNK/copy/skin/', recursive = F, full.names = F)

for(x in dirs){
  name <- x 
  #name <- gsub('_filtered_feature_bc_matrix','', x) 
  # x -> characters in dirs / _filtered_feature_bc_matrix -> expression pattern that the 
  # 'gsub' function is looking for in the string 'x'/ ' ' -> is the replacement string.  
  # if a pattern '_filtered_feature_bc_matrix' is found in 'x', it will be replaced with a space ' '
  
  
  cts <- ReadMtx(mtx = paste0('/home2/sbjs0428/R/atopyNK/copy/skin/',x,'/matrix.mtx.gz'),
                 features = paste0('/home2/sbjs0428/R/atopyNK/copy/skin/',x,'/features.tsv.gz'),
                 cells = paste0('/home2/sbjs0428/R/atopyNK/copy/skin/',x,'/barcodes.tsv.gz'))
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



obj <- merge(normalskin, y = c(atopyskin, atopyskintx),
                       add.cell.ids = c('NSK','ATSK', 'ATSKTX'),
                       project = 'ASK')

rm(dirs, atopyskin, atopyskintx, normalskin, name, cts)
# create a sample column
obj$sample <- rownames(obj@meta.data)

obj@meta.data <- separate(obj@meta.data, col = 'sample', into = c('Type', 'Barcode'), 
                                    sep = '_')



# calculate mitochondrial percentage
obj$mitoPercent <- PercentageFeatureSet(obj, pattern='^mt-')
obj <- subset(obj, subset = nCount_RNA > 500 &
                                   nFeature_RNA > 300 &
                                   mitoPercent < 10)

rm(sce,atopyskin, atopyskintx, normalskin, name, x)




# Calculate the number of cells in which each gene is expressed
gene_counts <- rowSums(obj[['RNA']]$counts > 0)

# Get genes that are expressed in at least 3 cells
genes_above_threshold <- names(gene_counts[gene_counts >= 3])


# Combine the genes above threshold with the genes in gene_list
final_gene_list <- genes_above_threshold

# Subset the Seurat object to retain only the genes in final_gene_list
obj <- subset(obj, features = final_gene_list)
rm(final_gene_list, genes_above_threshold, gene_counts)


# Extract all gene names
allgenes <- rownames(obj[['RNA']]$counts)

# Use grepl to create a logical vector for each condition
starts_with_Gm_or_Mt <- grepl("^(Gm|Mt|mt)", allgenes)
ends_with_Rik <- grepl("Rik$", allgenes)

# Combine conditions: genes that start with 'Gm' or 'Mt' OR end with 'Rik'
genes_to_exclude <- starts_with_Gm_or_Mt | ends_with_Rik

# Subset the allgenes vector to exclude these genes
filtered_genes <- allgenes[!genes_to_exclude]

allgenes_filtered <- filtered_genes
rm(starts_with_Gm_or_Mt, ends_with_Rik, genes_to_exclude, filtered_genes, allgenes)
obj <- subset(obj, features = allgenes_filtered)
rm(allgenes_filtered)



obj <- NormalizeData(obj)
obj <- FindVariableFeatures(object = obj)
obj <- ScaleData(object = obj, features = rownames(obj))
obj <- RunPCA(object = obj)
obj <- IntegrateLayers(obj, method =RPCAIntegration, orig.reduction='pca', new.reduction='pca', verbose=F)
obj <- FindNeighbors(object = obj, dims = 1:20, reduction='pca')
obj <- FindClusters(object = obj, resolution = 0.8, verbose = FALSE, cluster.name = 'seurat_clusters')
obj <- RunUMAP(object = obj, dims = 1:20, reduction = 'pca', reduction.name = 'umap')
obj <- JoinLayers(obj)


saveRDS(obj, '/home2/sbjs0428/R/atopyNK/workspace/yesbatch/skin_ref2018_obj.RDS')


