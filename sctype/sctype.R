library(Seurat)
library(tidyverse)
source('~/necessary_packages.R')
source('~/R/eQTU_rare_celltype/NARDsc/Rfunction/PracticalFunctionFindingRare.R')
rm(analyzerare, analyzetype, assignNAtoSmallClusters)
library(readxl)

#obj <- readRDS('../BPCells/cut_merged_NARD_res1.5_foundare.RDS')
#obj$sctype <- NA
# load libraries and functions
lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
source('/home2/sbjs0428/Tool/sc-type/R/mouse/gene_sets_prepare.R')

db_ = 'ScTypeDB_skin&lymph_common.xlsx'
#tissue = "PBMC" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 



#before run sctype function, object should be rescaled using all genes (In default, only 2000 variable features are scaled)


sctype_score <- function(scRNAseqData, scaled = !0, gs, gs2 = NULL, gene_names_to_uppercase = !0, ...){
  
  # check input matrix
  if(!is.matrix(scRNAseqData)){
    warning("scRNAseqData doesn't seem to be a matrix")
  } else {
    if(sum(dim(scRNAseqData))==0){
       warning("The dimension of input scRNAseqData matrix equals to 0, is it an empty matrix?")
    }
  }
  
  # marker sensitivity
  marker_stat = sort(table(unlist(gs)), decreasing = T); 
  marker_sensitivity = data.frame(score_marker_sensitivity = scales::rescale(as.numeric(marker_stat), to = c(0,1), from = c(length(gs),1)),
                                      gene_ = names(marker_stat), stringsAsFactors = !1)

#  # convert gene names to Uppercase
#  if(gene_names_to_uppercase){
#    rownames(scRNAseqData) = toupper(rownames(scRNAseqData));
#  }
  
  # subselect genes only found in data
  names_gs_cp = names(gs); names_gs_2_cp = names(gs2);
  gs = lapply(1:length(gs), function(d_){ 
    GeneIndToKeep = rownames(scRNAseqData) %in% as.character(gs[[d_]]); rownames(scRNAseqData)[GeneIndToKeep]})
  gs2 = lapply(1:length(gs2), function(d_){ 
    GeneIndToKeep = rownames(scRNAseqData) %in% as.character(gs2[[d_]]); rownames(scRNAseqData)[GeneIndToKeep]})
  names(gs) = names_gs_cp; names(gs2) = names_gs_2_cp;
  cell_markers_genes_score = marker_sensitivity[marker_sensitivity$gene_ %in% unique(unlist(gs)),]
  
  # z-scale if not
  if(!scaled) Z <- t(scale(t(scRNAseqData))) else Z <- scRNAseqData
  
  # multiple by marker sensitivity
  for(jj in 1:nrow(cell_markers_genes_score)){
    Z[cell_markers_genes_score[jj,"gene_"], ] = Z[cell_markers_genes_score[jj,"gene_"], ] * cell_markers_genes_score[jj, "score_marker_sensitivity"]
  }
  
  # subselect only with marker genes
  Z = Z[unique(c(unlist(gs),unlist(gs2))), ]
  
  # combine scores
  es = do.call("rbind", lapply(names(gs), function(gss_){ 
    sapply(1:ncol(Z), function(j) {
      gs_z = Z[gs[[gss_]], j]; gz_2 = Z[gs2[[gss_]], j] * -1
      sum_t1 = (sum(gs_z) / sqrt(length(gs_z))); sum_t2 = sum(gz_2) / sqrt(length(gz_2));
      if(is.na(sum_t2)){
        sum_t2 = 0;
      }
      sum_t1 + sum_t2
    })
  })) 
  
  dimnames(es) = list(names(gs), colnames(Z))
  es.max <- es[!apply(is.na(es) | es == "", 1, all),] # remove na rows
 
  es.max
}

sctype <- function(obj, tissue) {
        gs_list = gene_sets_prepare(db_, tissue)
         # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 
        if (tissue == 'auto') {
                source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/auto_detect_tissue_type.R")
                tissue = auto_detect_tissue_type(path_to_db_file = db_, seuratObject = obj, scaled = TRUE, assay = "RNA")
        } else {
        }
        # prepare gene sets
        gs_list = gene_sets_prepare(db_, tissue)

        es.max = sctype_score(scRNAseqData = obj[["RNA"]]$scale.data, scaled = T,
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
        cL_resutls = do.call("rbind", lapply(unique(obj@meta.data$seurat_clusters), function(cl){
        es.max.cl = sort(rowSums(es.max[ ,rownames(obj@meta.data[obj@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
        head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(obj@meta.data$seurat_clusters==cl)), 10)
}))
        sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

        # set low-confident (low ScType score) clusters to "unknown"
        sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
        obj@meta.data$sctype = ''
        for(j in unique(sctype_scores$cluster)){
                cl_type = sctype_scores[sctype_scores$cluster==j,];
                obj@meta.data$sctype[obj@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
        }
        return(obj)
}


for (i in c('CD4 T cell','CD8 T cell', 'Dendritic cell', 'Granulocytes','B cell', 'NK')) {
        if (i %in% unique(diff$sctype)) {
                temp <- subset(diff, sctype == i)
                temp <- NormalizeData(temp)
                temp <- FindVariableFeatures(temp)
                temp <- ScaleData(temp, features = rownames(temp))
                temp <- sctype(temp, tissue=i)
                if ('Unknown' %in% unique(temp$sctype)) {
                        temp$sctype[temp$sctype == 'Unknown'] <- paste0(i, ' ', 'notclassified')
                diff@meta.data <- diff@meta.data %>%
                left_join(temp@meta.data %>% select(barcode, sctype), by = "barcode", suffix = c("", ".temp")) %>%
                mutate(sctype = coalesce(sctype.temp, sctype)) %>%
                select(-sctype.temp)  # Clean up the temporary column
                rownames(diff@meta.data) <- diff$barcode
                }
        }
}

