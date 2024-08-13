source('~/necessary_packages.R')

library(gage)
library(pathview)
library(gageData)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(rWikiPathways)
library(KEGGREST)





load('/gpfs/home2/sbjs0428/R/atopyNK/workspace/nobatch/lymph_latest.RData')


path_list <- keggList('pathway', 'mmu')

obj_list <- list()

deg_df <- data.frame(row.names = levels(obj$celltype), CellType = levels(obj$celltype), NumDEGs = 0)
kegg_df <- data.frame(row.names = levels(obj$celltype), CellType = levels(obj$celltype), NumKEGGs = 0, NumEnKEGGs = 0)
wiki_df <- data.frame(row.names = levels(obj$celltype), CellType = levels(obj$celltype), NumWKs = 0)


for (ct in levels(obj$celltype)) {
        obj_list[[ct]] <- list()
        temp_obj <- subset(obj, celltype == ct)
        obj_list[[ct]][['obj']] <- temp_obj
        tryCatch({
                temp_deg <- FindMarkers(temp_obj, group.by='Type', ident.1='ATLNTX', ident.2='ATLN') %>% filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.4)
                temp_deg$gene <- rownames(temp_deg)
                temp_id <- bitr(rownames(temp_deg), fromType='SYMBOL', toType='ENTREZID', OrgDb=org.Mm.eg.db)
                gene_indices <- match(temp_deg$gene, temp_id$SYMBOL)
                temp_deg$id <- temp_id[gene_indices, 'ENTREZID']
                rm(temp_id, gene_indices)
                library(openxlsx)
                write.xlsx(temp_deg, paste0('/gpfs/home2/sbjs0428/R/atopyNK/deg_xlsx/nobatch/lymph/', ct, '_deg.xlsx'), asTable=T)
                obj_list[[ct]][['deg']] <- temp_deg
                deg_df[ct, 'NumDEGs'] <- nrow(temp_deg)
                obj_list[[ct]][['input']] <- kegginput(temp_deg)
        }, error = function(e) {
                message('There are no DEGs for cell type: ', ct)
                temp_deg <<- 'No significant DEG'
                obj_list[[ct]][['deg']] <<- temp_deg
        })
        tryCatch({
                if (!is.character(temp_deg)) {
                        enrich_kegg <- enrichKEGG(names(kegginput(temp_deg)), organism='mmu', use_internal_data=T)@result
                        enrich_kegg <- enrich_kegg %>% arrange(desc(Count))
                        enrich_kegg$Description <- sapply(enrich_kegg$ID, function(id) path_list[[id]])
                        enrich_kegg <- enrich_kegg %>% filter(Count >= 5)
                        kegg_df[ct, 'NumEnKEGGs'] <- nrow(enrich_kegg)
                        if (nrow(enrich_kegg) > 0) {
                                obj_list[[ct]][['enrichkegg']] <- enrich_kegg
                                print(paste('Enriched KEGG pathways for cell type', ct, ':', nrow(enrich_kegg)))
                                write.xlsx(enrich_kegg, paste0('/home2/sbjs0428/R/atopyNK/kegg/lymph/kegg_xlsx/', ct, '_enrichedKEGG.xlsx'), asTable=T)      
                        } else {
                                obj_list[[ct]][['enrichkegg']] <- 'No significant enrich KEGG'
                        }
                }
        }, error = function(e) {
                message('There are no significant enrich KEGG for cell type: ', ct)
                enrich_kegg <<- 'No significant enrich KEGG'
                obj_list[[ct]][['enrichkegg']] <<- enrich_kegg
        })
        tryCatch({
                if (!is.character(temp_deg)) {
                        temp_kegg <- deg2kegg(kegginput(temp_deg), celltype = ct)
                        if (nrow(temp_kegg) > 0) {
                                obj_list[[ct]][['kegg']] <- temp_kegg
                                kegg_df[ct, 'NumKEGGs'] <- nrow(temp_kegg)
                                print(paste('KEGG pathways for cell type', ct, ':', nrow(temp_kegg)))
                        } else {
                                obj_list[[ct]][['kegg']] <- 'No significant KEGG'
                        }
                } else {
                        obj_list[[ct]][['kegg']] <- 'No significant KEGG'
                        message('There are no KEGG as no DEGs were found')
                }
        }, error = function(e) {
                message('There are no significant KEGG for cell type: ', ct)
                temp_kegg <<- 'No significant KEGG'
                obj_list[[ct]][['kegg']] <<- temp_kegg
        })
        tryCatch({
                if (!is.character(temp_deg)) {
                        temp_wiki <- enrichWP(temp_deg$id, organism='Mus musculus')@result
                        temp_wiki <- temp_wiki  %>% arrange(desc(Count))
                        if (nrow(temp_wiki) > 0) {
				wiki_df[ct, 'NumWKs'] <- nrow(temp_wiki)
                                obj_list[[ct]][['wiki']] <- temp_wiki
                                print(paste('WIKI pathways for cell type', ct, ':', nrow(temp_wiki)))
                                write.xlsx(temp_wiki, paste0('/home2/sbjs0428/R/atopyNK/wikipath/lymph/wiki_xlsx/', ct, '_wiki.xlsx'), asTable=T)
                                cat('\n\n')
                        } else {
                                obj_list[[ct]][['wiki']] <- 'No significant WIKI'
                                cat('\n\n')
                        }
                } else {
                        obj_list[[ct]][['wiki']] <- 'No significant WIKI'
                        cat('\n\n')
                }
        }, error = function(e) {
                message('There are no significant WIKI for cell type: ', ct)
                obj_list[[ct]][['wiki']] <<- 'No significant WIKI'
                cat('\n\n')
                })
        if (ct == levels(obj$celltype)[[length(levels(obj$celltype))]]) {
                rm(ct, temp_kegg, temp_obj, temp_deg, enrich_kegg, temp_wiki)
        }
}


save.image('/gpfs/home2/sbjs0428/R/atopyNK/workspace/nobatch/lymph_latest.RData')
