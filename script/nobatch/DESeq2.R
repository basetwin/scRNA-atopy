deseq <- function(obj, Type = 'Type', ref, alpha) {
    # Extract counts and metadata
    cts <- obj[['RNA']]$counts
    colData <- obj@meta.data
    Type_formula <- as.formula(paste("~ ", Type))
    # Create a DESeqDataSet
    dds <- DESeqDataSetFromMatrix(countData = cts, colData = colData, design = Type_formula)
    keep <- rowSums(counts(dds)) >= 1  
    dds <- dds[[keep,]	
    dds[[Type]] <- relevel(dds[[Type]], ref = ref)
    dds <- DESeq(dds)
    
    if (length(alpha) > 1) {
        # Initialize a list to store results if multiple alpha values are provided
        res_list <- list()
        for (num in alpha) {
            res_name <- paste0("res", gsub("\\.", "", format(num)))
            res_list[[res_name]] <- results(dds, alpha=num)
        }
        return(list(result = res_list, dds = dds))
    } else {
        return(list(result = results(dds, alpha=alpha), dds = dds))
    }
}

#        res_name <- paste0("res", gsub("\\.", "", format(alpha)))
#        res <- results(dds, alpha = alpha)
#        assign(res_name, res, envir = .GlobalEnv)
#res0.01 <- results(dds, alpha = 0.01)
#summary(res0.01)
#
## contrasts
#resultsNames(dds)
#
## e.g.: treated_4hrs, treated_8hrs, untreated
#
#results(dds, contrast = c("dexamethasone", "treated_4hrs", "untreated"))
#
## MA plot
#plotMA(res)
