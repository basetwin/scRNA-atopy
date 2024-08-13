deg_list <- list()  # Initialize the list to store the results

celltypes <- levels(obj$celltype)

for (i in celltypes) {
    tryCatch({
        # Subset the object to the specific cell type
        subset_obj <- subset(obj, subset = celltype == i)

        # Run FindMarkers and filter the results
        markers <- FindMarkers(subset_obj, group.by = 'Type', ident.1 = 'ATSKTX', ident.2 = 'ATSK') %>%
                   filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.6)

        # Get the number of rows
        num <- nrow(markers)

        # Assign the result to the list
        deg_list[paste0('Celltype: ', i)] <- num
    }, error = function(e) {
        # In case of an error, assign 0
        deg_list[paste0('Celltype: ', i)] <- 0
        message(paste("Error for cell type", i, ":", e$message))  # Print error message for debugging
    })}


