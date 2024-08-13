z_score_celltype <- function(obj, genes) {
  library(dplyr)
  library(Seurat)
  gene_exp <- FetchData(obj, vars = genes, slot = 'counts')
  meta <- obj@meta.data
  
  z_scores_list <- list()
  
  for (gene in genes) {
    gene_exp_vector <- gene_exp[[gene]]
    
    mean_g <- mean(gene_exp_vector, na.rm = TRUE)
    sd_g <- sd(gene_exp_vector, na.rm = TRUE)
    
    z_scores <- (gene_exp_vector - mean_g) / sd_g
    meta[[paste0("z_score_", gene)]] <- z_scores
    
    z_score_clusters <- meta %>%
      group_by(seurat_clusters) %>%
      summarize(mean_z_score = mean(get(paste0("z_score_", gene)), na.rm = TRUE))
    
    z_scores_list[[gene]] <- z_score_clusters
  }
  z_scores_df <- Reduce(function(x, y) full_join(x, y, by = "seurat_clusters"), z_scores_list)
  z_scores_df$seurat_clusters <- NULL
  colnames(z_scores_df) <- genes 
  z_scores_mean <- data.frame(mean_zscore = rowMeans(z_scores_df)) 
  return(list(z_scores_pergene = z_scores_df, z_score_mean = z_scores_mean))
}
