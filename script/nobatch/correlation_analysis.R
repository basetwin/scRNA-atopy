# Load the Seurat object
obj <- readRDS('/gpfs/home2/sbjs0428/R/atopyNK/workspace/nobatch/skin_obj_240716.RDS')

# Aggregate expression data
pseudo_obj <- AggregateExpression(obj, return.seurat = TRUE, assays = 'RNA', group.by = c('Type', 'celltype'))

# Split the Seurat object by cell type
obj_list <- SplitObject(pseudo_obj, split.by = 'celltype')

# Define the function to find DEGs and calculate correlation matrices
findeg <- function(obj, condition1, condition2) {
  deg <- FindMarkers(obj, group.by = 'Type', ident.1 = condition1, ident.2 = condition2) %>% 
    filter(p_val_adj < 0.05)
  return(rownames(deg))}

calc_correlation <- function(pseudo_obj, genes, condition1, condition2) {
  data1 <- as.data.frame(pseudo_obj@assays$RNA$counts[genes, pseudo_obj$Type == condition1])
  data2 <- as.data.frame(pseudo_obj@assays$RNA$counts[genes, pseudo_obj$Type == condition2])
  cor_matrix <- cor(data1, data2)
  return(cor_matrix)}

# Initialize the results list
cor_results <- list()

# Calculate correlations for each cell type
for (celltype in names(obj_list)) {
  object <- subset(obj, celltype == celltype)
  genes <- findeg(object, 'NSK', 'ATSK')
  pseudo_subset <- obj_list[[celltype]]
  cor_norm_atopic <- calc_correlation(pseudo_subset, genes, "NSK", "ATSK")
  cor_norm_treated <- calc_correlation(pseudo_subset, genes, "NSK", "ATSKTX")
  cor_results[[celltype]] <- list(norm_atopic = cor_norm_atopic, norm_treated = cor_norm_treated)
}

# Function to plot correlation matrices
plot_correlations <- function(cor_results, celltype) {
  cor_norm_atopic <- cor_results[[celltype]]$norm_atopic
  cor_norm_treated <- cor_results[[celltype]]$norm_treated
  df_norm_atopic <- melt(cor_norm_atopic)
  df_norm_treated <- melt(cor_norm_treated)
  
  p1 <- ggplot(df_norm_atopic, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
    theme_minimal() +
    ggtitle(paste("Normal vs Atopic Skin -", celltype))
  
  p2 <- ggplot(df_norm_treated, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
    theme_minimal() +
    ggtitle(paste("Normal vs Treated Atopic Skin -", celltype))
  
  grid.arrange(p1, p2, nrow = 1)}

# Plot correlations for each cell type
for (celltype in names(obj_list)) {
  plot_correlations(cor_results, celltype)}
















######################################################
######################################################



# Load the Seurat object
obj <- readRDS('/gpfs/home2/sbjs0428/R/atopyNK/workspace/nobatch/skin_obj_240716.RDS')

# Aggregate expression data
pseudo_obj <- AggregateExpression(obj, return.seurat = TRUE, assays = 'RNA', group.by = c('Type', 'celltype'))

# Split the Seurat object by cell type
obj_list <- SplitObject(pseudo_obj, split.by = 'celltype')

# Define the function to find DEGs
findeg <- function(obj, condition1, condition2) {
  deg <- FindMarkers(obj, group.by = 'Type', ident.1 = condition1, ident.2 = condition2) %>%       
    filter(p_val_adj < 0.05)
  return(rownames(deg))
}

# Initialize the results list
cor_results <- list()

# Calculate correlations and prepare data for each cell type
for (celltype in names(obj_list)) {
  object <- subset(obj, celltype == celltype)
  genes <- findeg(object, 'NSK', 'ATSK')
  pseudo_subset <- obj_list[[celltype]]
  data1 <- as.data.frame(pseudo_subset@assays$RNA$counts[genes, pseudo_subset$Type == "NSK"])
  data2 <- as.data.frame(pseudo_subset@assays$RNA$counts[genes, pseudo_subset$Type == "ATSK"])
  data3 <- as.data.frame(pseudo_subset@assays$RNA$counts[genes, pseudo_subset$Type == "ATSKTX"])
  data_combined <- list(NSK = data1, ATSK = data2, ATSKTX = data3)
  cor_results[[celltype]] <- data_combined
}


plot_scatter_matrix <- function(data, celltype) {
  data_combined <- do.call(cbind, data)
  condition_labels <- rep(names(data), sapply(data, ncol))
  colnames(data_combined) <- paste(condition_labels, 1:ncol(data_combined), sep = "_")
  
  if (length(data) == 3) {
    ggpairs(data_combined,
            upper = list(continuous = wrap("cor", size = 5)),
            lower = list(continuous = wrap("points", alpha = 0.3, size = 0.5))) +
      ggtitle(paste("Correlation between NSK, ATSK, and ATSKTX -", celltype))
  } else if (length(data) == 2) {
    pairs(data_combined,
          main = paste("Pairwise Plot -", celltype),
          pch = 19, col = alpha("blue", 0.3))
  } else {
    ggplot(melt(data_combined), aes(x = Var2, y = value)) +
      geom_point(alpha = 0.3) +
      ggtitle(paste("Single Condition Plot -", celltype))
  }}



# Plot scatter plot matrix for each cell type
for (celltype in names(cor_results)) {
  plot_scatter_matrix(cor_results[[celltype]], celltype)
}
