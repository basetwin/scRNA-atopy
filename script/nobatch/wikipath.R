

##  pathway archive should be downloaded passively at http://data.wikipathways.org -> current -> mus musculus 


wp2gene <- as.data.frame(readPathwayGMT('/home2/sbjs0428/R/atopyNK/wikipath/wikipathways-20240710-gmt-Mus_musculus.gmt'))


wpid2gene <- wp2gene %>% dplyr::select(wpid, gene)
wpid2name <- wp2gene %>% dplyr::select(wpid, name)


#all_ids <- deg$id
#up_ids <- deg[deg$avg_log2FC > 0.4, 'id']
#dn_ids <- deg[deg$avg_log2FC < -0.4, 'id']
#ewp <- enricher( up_ids, universe=deg$id, OrgDb = org.Mm.eg.db, readable=T)



wikienrich <- function(celltype) {
	deg <- obj_list[[celltype]]$deg
	all_ids <- deg$id
	up_ids <- deg[deg$avg_log2FC > 0.4, 'id']
	dn_ids <-deg[deg$avg_log2FC < -0.4, 'id']
	ewp <- enricher(gene = up_ids, universe = all_ids, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
	ewp <- DOSE::setReadable(ewp, org.	Mm.eg.db, keyType='ENTREZID')
	#ewp <- ewp@result %>% arrange(desc(Count))
	return(ewp)
}


wikifull <- function(celltype) {
	deg <- obj_list[[celltype]]$deg
	ewp <- enrichWP(deg$id, organism= 'Mus musculus')
	ewp <- ewp@result
	
