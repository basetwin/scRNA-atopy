
library(KEGGREST)
library(gage)
library(pathview)
library(gageData)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(rWikiPathways)



## this does not work for our data as our data has gene official name, but pathview need entrez id of gene so we should change our gene name into entrez id using the code below

deg2kegg <- function(kegg_gene_list, organism='mouse', celltype = 'example') { 
	library(gage)
	library(pathview)
	library(gageData)
	if (organism == 'mouse') {
		library(org.Mm.eg.db)
		database = org.Mm.eg.db
		species = 'mmu'
	} else if (organism == 'human') {
		library(org.Hs.eg.db)
		database = org.Hs.eg.db
		species = 'hsa'
	}
	if (celltype == 'example') {
		stop('Please assign celltype')
	}
	path_list <- keggList('pathway', species)
	library(KEGG.db)
	kk <- gseKEGG(geneList = kegg_gene_list, organism= species, keyType ='kegg', use_internal_data=T, pvalueCutoff=1)
	enriched_path <- path_list[ kk@result$ID]
	enriched_path <- as.data.frame(enriched_path)
	enriched_path$pvalue <- kk@result$pvalue
	enriched_path$id <- rownames(enriched_path) 
	return(enriched_path)
}


kegginput <- function(deg, organism='mouse') {
	if (organism == 'mouse') {
		library(org.Mm.eg.db)
		database = org.Mm.eg.db
                species = 'mmu'
	} else if (organism == 'human') {
                library(org.Hs.eg.db)
                database = org.Hs.eg.db
                species = 'hsa'
	}
	kegg_gene_list <- deg$avg_log2FC
        library(KEGGREST)
	genes <- mapIds(database,keys = rownames(deg), column = 'ENTREZID', keytype = 'SYMBOL', multiVals='first')
	names(kegg_gene_list) <- genes
	kegg_gene_list <- kegg_gene_list[ !is.na(names(kegg_gene_list))]
	kegg_gene_list <- sort(kegg_gene_list, decreasing=T)
	return(kegg_gene_list)
}

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
		temp_deg <- FindMarkers(temp_obj, group.by='Type', ident.1='ATSKTX', ident.2='ATSK') %>% filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.4)
		temp_deg$gene <- rownames(temp_deg)
		temp_id <- bitr(rownames(temp_deg), fromType='SYMBOL', toType='ENTREZID', OrgDb=org.Mm.eg.db)
		gene_indices <- match(temp_deg$gene, temp_id$SYMBOL)
		temp_deg$id <- temp_id[gene_indices, 'ENTREZID']
		rm(temp_id, gene_indices)
		library(openxlsx)
		write.xlsx(temp_deg, paste0('/gpfs/home2/sbjs0428/R/atopyNK/deg_xlsx/nobatch/skin/', ct, '_deg.xlsx'), asTable=T)
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
				write.xlsx(enrich_kegg, paste0('/home2/sbjs0428/R/atopyNK/kegg/skin/kegg_xlsx/', ct, '_enrichedKEGG.xlsx'), asTable=T)	
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
				write.xlsx(temp_wiki, paste0('/home2/sbjs0428/R/atopyNK/wikipath/skin/wiki_xlsx/', ct, '_wiki.xlsx'), asTable=T)
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


############# using biocarta enriching 

biocarta <- function(deg) {
	







for (g in rownames(cts)) {
min_val <- min(cts[g, ])
max_val <- max(cts[g, ])
val <- c() 
for (sam in colnames(cts)) {
val <- c(val, (cts[g, sam] - min_val) / (max_val - min_val))
}
cts[g, ] <- val
}






# a <- deg_df$CellType[ deg_df$NumDEGs == 0]
#for (i in a) {
#
#tem <- subset(obj, celltype == i)
#
#tryCatch({
#
#a <- FindMarkers(tem, group.by='Type', ident.1='ATSKTX', ident.2='ATSK') %>% filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.4)
#print(paste0(i, ' type deg num: ', nrow(a)))
#print(rownames(a))
#}, error = function(e) {
#
#message('No degs')
#
#})}

'''

"#FC8D59" "#FFFFBF" "#91BFDB"
   red       yellow     blue


'''
savepath <- function(celltype, pathid, type='skin') {
	current_path <- getwd()
	geneList <- obj_list[[celltype]]$input
	dir_name <- gsub(' ', '', celltype)
	dir_name <- paste0('/gpfs/home2/sbjs0428/R/atopyNK/kegg/', type, '/keggimage/', dir_name)
	if (!dir.exists(dir_name)) {
		dir.create(dir_name, recursive=T)
	}

	setwd(dir_name)
	pathview(geneList, species='mmu', pathway.id = pathid, limit = list(gene=max(abs(geneList))), low=list(gene='#91BFDB'), mid=list(gene='#FFFFBF'), high=list(gene='#FC8D59'))
	files <- list.files(dir_name, full.names=T)
	files_to_remove <- files[!grepl('pathview', files)]
	if (length(files_to_remove) > 0) {
		file.remove(files_to_remove)
	}
	setwd(current_path)
}


ggplot(deg_df, aes(x = reorder(CellType, NumDEGs), y = NumDEGs)) +
  geom_bar(stat = "identity", alpha=0.6, aes(fill=CellType)) +
  geom_text(aes(label = NumDEGs),position = position_stack(vjust = 0.5),hjust = -0.2) +
  coord_flip() +
  labs(title = "Number of DEGs across Cell Types (FC > 1.3)", x = "Cell Type", y = "Number of DEGs") +
  theme_minimal()



ggplot(kegg_df, aes(x= reorder(CellType, NumKEGGs), y=NumKEGGs)) +
	geom_bar(stat = 'identity', aes(fill=CellType), alpha=0.6) +
	geom_text(aes(label = NumKEGGs),position = position_stack(vjust = 0.5),hjust = -0.2) +
	coord_flip() +
	labs(title = 'Number of KEGGs across Cell Types', x='Cell Type', y='Number of KEGGs') +
	theme_minimal()

ggplot(kegg_df, aes(x= reorder(CellType, NumEnKEGGs), y=NumEnKEGGs)) +
  labs(title = "Number of DEGs across Cell Types (FC > 1.3)", x = "Cell Type", y = "Number of DEGs") +
  theme_minimal()



ggplot(kegg_df, aes(x= reorder(CellType, NumKEGGs), y=NumKEGGs)) +
	geom_bar(stat = 'identity', aes(fill=CellType), alpha=0.6) +
	geom_text(aes(label = NumKEGGs),position = position_stack(vjust = 0.5),hjust = -0.2) +
	coord_flip() +
	labs(title = 'Number of KEGGs across Cell Types', x='Cell Type', y='Number of KEGGs') +
	theme_minimal()

ggplot(kegg_df, aes(x= reorder(CellType, NumEnKEGGs), y=NumEnKEGGs)) +
        geom_bar(stat = 'identity', aes(fill=CellType), alpha=0.6) +
	geom_text(aes(label = NumEnKEGGs),position = position_stack(vjust = 0.5),hjust = -0.2) +
        coord_flip() +
        labs(title = 'Number of enriched KEGGs across Cell Types', x='Cell Type', y='Number of KEGGs') +
        theme_minimal()

pathplot <- function(path) {
	if (! path %in% c('deg', 'kegg', 'wiki')) {
		stop('Variable path should be one of "deg", "kegg", "wiki"')
	}
	if (path == 'wiki') {
		ggplot(wiki_df, aes(x= reorder(CellType, NumWKs), y=NumWKs)) +
		geom_bar(stat = 'identity', aes(fill=CellType), alpha=0.6) +
		geom_text(aes(label = NumWKs),position = position_stack(vjust = 0.5),hjust = -0.2) +
		coord_flip() +
		labs(title = 'Number of enriched WIKIs across Cell Types', x='Cell Type', y='Number of WIKIs') +
		theme_minimal()
	} else if (path == 'kegg') {
		ggplot(kegg_df, aes(x= reorder(CellType, NumEnKEGGs), y=NumEnKEGGs)) +
		geom_bar(stat = 'identity', aes(fill=CellType), alpha=0.6) +
		geom_text(aes(label = NumEnKEGGs),position = position_stack(vjust = 0.5),hjust = -0.2) +
		coord_flip() +
		labs(title = 'Number of enriched KEGGs across Cell Types', x='Cell Type', y='Number of KEGGs') +
		theme_minimal()
	} else {
		ggplot(deg_df, aes(x = reorder(CellType, NumDEGs), y = NumDEGs)) +
		geom_bar(stat = "identity", alpha=0.6, aes(fill=CellType)) +
		geom_text(aes(label = NumDEGs),position = position_stack(vjust = 0.5),hjust = -0.2) +
		coord_flip() +
		labs(title = "Number of DEGs across Cell Types (FC > 1.3)", x = "Cell Type", y = "Number of DEGs") +
		theme_minimal()
	}
}

for (i in c('deg', 'kegg', 'wiki')) {
	plot <- pathplot(i)
	ggsave(paste0('/gpfs/home2/sbjs0428/R/atopyNK/summaryplot/skin/', i, '_summary.png'), plot, width=10, height=6, units='in', dpi=1200, bg='white')
}

ggplot(wiki_df, aes(x= reorder(CellType, NumWKs), y=NumWKs)) +
        geom_bar(stat = 'identity', aes(fill=CellType), alpha=0.6) +
	geom_text(aes(label = NumWKs),position = position_stack(vjust = 0.5),hjust = -0.2) +
        coord_flip() +
        labs(title = 'Number of enriched WIKIs across Cell Types', x='Cell Type', y='Number of WIKIs') +
        theme_minimal()

################################################################################


fibro <- subset(obj, celltype =='Fibroblast')

fibro_paths <- c('04657', '04668', '04060','04062', '04064', '04659')

deg <- FindMarkers(fibro, group.by='Type', ident.1='ATSKTX', ident.2 = 'ATSK') %>% filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.6)

kegg_gene_list <- kegginput(deg)

celltype = 'fibroblast'

for (path in fibro_paths) {
	if (! dir.exists(paste0('/home2/sbjs0428/R/atopyNK/kegg/', celltype, '/'))) {
		dir.create(paste0('/home2/sbjs0428/R/atopyNK/kegg/', celltype, '/'))
	}
	pathview(kegg_gene_list, pathway.id = paste0('mmu', path), species='mmu', kegg.dir = paste0('/home2/sbjs0428/R/atopyNK/kegg/', celltype, '/'), res=1200)
}	

gk <- subset(obj, celltype = 'G keratinocyte')

gk_paths <- c('04657', '04060', '04660', '04061', '04650', '04151', '05235', '04630', '04750', '04010', '04020', '04068', '04659', '04668', '04062', '04530', '04064', '04510', '04072', '04670', '04919', '04310', '04620', '04066')

findeg <- function(obj) { 
deg <- FindMarkers(obj, group.by='Type', ident.1='ATSKTX', ident.2='ATSK') %>% filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.6)
return(deg)
}

deg <-  findeg(gk)

kegg_gene_list <- kegginput(deg)

celltype = 'gkera'

for (path in gk_paths) {
        if (! dir.exists(paste0('/home2/sbjs0428/R/atopyNK/kegg/', celltype, '/'))) {
                dir.create(paste0('/home2/sbjs0428/R/atopyNK/kegg/', celltype, '/'))
        }
        pathview(kegg_gene_list, pathway.id = paste0('mmu', path), species='mmu', kegg.dir = paste0('/home2/sbjs0428/R/atopyNK/kegg/', celltype, '/'), res=300)
}


a <- c('mmu04062', 'mmu04310', 'mmu04660', 'mmu04330', 'mmu04630', 'mmu04350', 'mmu04670', 'mmu04664', 'mmu04659', 'mmu04658', 'mmu04668', 'mmu01521', 'mmu04657', 'mmu04060', 'mmu04064')
for (path in a) {
	pathview(foldchanges, pathway.id = path, species = 'mmu', kegg.dir  = '/home2/sbjs0428/R/atopyNK/kegg')
} 



#mmu04062 Chemokine signaling pathway
#mmu04310 Wnt signaling pathway
#mmu04660 T cell receptor signaling pathway
#mmu04330 Notch signaling pathway
#mmu04630 Jak-STAT signaling pathway
#mmu04350 TGF-beta signaling pathway
#mmu04670 Leukocyte transendothelial migration
#mmu04664 Fc epsilon RI signaling pathway
#mmu04659 Th17 cell differentiation - Mus musculus (house mouse)
#mmu04658 Th1 and Th2 cell differentiation - Mus musculus (house mouse)
#mmu04668 TNF signaling pathway - Mus musculus (house mouse)
#mmu01521 EGFR tyrosine kinase inhibitor resistance - Mus musculus (house mouse)
#mmu04657 IL-17 signaling pathway - Mus musculus (house mouse)
#mmu04060 Cytokine-cytokine receptor interaction - Mus musculus (house mouse)
#mmu04064 NF-kappa B signaling pathway - Mus musculus (house mouse)




data('go.sets.mm')
data('go.subs.mm')

gobpsets <- go.sets.mm[go.subs.mm$BP]
gobpres <- gage(exprs = foldchanges, gsets = gobpsets, same.dir=T)

data(kegg.sets.mm)
data(sigmet.idx.mm)

kegg.sets.mm <- kegg.sets.mm[sigmet.idx.mm]

kegres <- gage(exprs = foldchanges, gsets = kegg.sets.mm, same.dir=T)

keggrespathways <- data.frame(id = rownames(kegres$greater), kegres$greater) %>%
					tibble::as_tibble() %>% 
					filter(row_number() <= 20) %>%
					.$id %>% as.character()

keggresids = substr(keggrespathways, start=1, stop=8)

tmp = sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species = 'mmu'))










## making local database

remotes::install_github("YuLab-SMU/createKEGGdb")
## clusterProfiler should be detached 

library(createKEGGdb)
species <-c("ath","hsa","mmu", "rno","dre","dme","cel")
createKEGGdb::create_kegg_db(species)
install.packages("KEGG.db_1.0.tar.gz", repos=NULL,type="source")
library(KEGG.db)

## up code works for me



# Produce the native KEGG plot (PNG)
dme <- pathview(gene.data=kegg_gene_list, pathway.id="dme04130", species = kegg_organism)

# Produce a different plot (PDF) (not displayed here)
dme <- pathview(gene.data=kegg_gene_list, pathway.id="dme04130", species = kegg_organism, kegg.native = F)




################################################################################################################



### for Reactome pathway analysis

library(ReactomePA)

deg_id <- kegginput(deg)

ra_list <- gsePathway(deg_id, pvalueCutoff=1, organism='mouse')@result





######## for wikipath analysis

library(RCy3)

cytoscapePing()

installApp('WikiPathways')
commandsRun(paste0('wikipathways import-as-pathway id=', pathid)) 

##  pathway archive should be downloaded passively at http://data.wikipathways.org -> current -> mus musculus 


wp2gene <- as.data.frame(readPathwayGMT('/home2/sbjs0428/R/atopyNK/wikipath/wikipathways-20240710-gmt-Mus_musculus.gmt'))
wikienrich <- function(celltype) {
        deg <- obj_list[[celltype]]$deg
        all_ids <- deg$id
        up_ids <- deg[deg$avg_log2FC > 0.4, 'id']
        dn_ids <-deg[deg$avg_log2FC < -0.4, 'id']
        ewp <- enricher(gene = up_ids, universe = all_ids, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
        ewp <- DOSE::setReadable(ewp, org.      Mm.eg.db, keyType='ENTREZID')
        #ewp <- ewp@result %>% arrange(desc(Count))
        return(ewp)
}


wikifull <- function(celltype) {
        deg <- obj_list[[celltype]]$deg
        ewp <- enrichWP(deg$id, organism= 'Mus musculus')
        ewp <- ewp@result
	return(ewp)
}
