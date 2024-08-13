source('~/necessary_packages.R')


subsampling <- function(tissue = 'skin') {
	if (!tissue %in% c('skin', 'lymph')) {
		stop("Please set the tissue type to 'skin' or 'lymph'")
	}	
	if (tissue == 'skin') {
		print('Tissue type is skin')
		obj <- readRDS('/home2/sbjs0428/R/atopyNK/workspace/nobatch/skin_obj_240716.RDS')
		stand_cell <- 'Bulge cell'
		objlist_type <- SplitObject(obj, split.by='Type')
		st_num <- unlist(lapply(objlist_type, function(x) table(x@meta.data$celltype)[stand_cell]))
		st_num_min <- min(st_num)
		print(sprintf('Mininum number of standard cell type is %d', st_num_min))
	} else {
		print('Tissue type is lymph')
		obj <- readRDS('/home2/sbjs0428/R/atopyNK/workspace/nobatch/lymph_obj.RDS')
		objlist_type <- SplitObject(obj, split.by='Type')
		stand_cell <- 'Naive CD8 T cell'
		st_num_min <- lapply(objlist_type, function(x) table(x@meta.data$celltype)[[stand_cell]]) %>% unlist() %>% min()
		print(sprintf('Mininum number of standard cell type is %d', st_num_min))
	}
	sampling_names <- c()
	for (type in names(objlist_type)) {
		tmp_obj <- objlist_type[[type]]
		sampling_ratio <- st_num_min / table(tmp_obj@meta.data$celltype)[[stand_cell]]
		print(sprintf('Sampling ration of %s is %f', type, sampling_ratio)) 
		type_bycell <- SplitObject(tmp_obj, split.by = 'celltype')
		for (cell in names(type_bycell)) {
			if (cell == stand_cell) {
				cell_names <- colnames(type_bycell[[cell]])
				cell_names <- sample(cell_names, size=st_num_min)
				print(sprintf("%s's number is %d", cell, length(cell_names)))
				sampling_names <- c(cell_names, sampling_names)
				#type_bycell[[cell]] <- subset(type_bycell[[cell]], cells = cell_names)
			} else {
				sampling_num <- ncol(type_bycell[[cell]]) * sampling_ratio
				cell_names <- colnames(type_bycell[[cell]])
				cell_names <- sample(cell_names, size = sampling_num)
				sampling_names <- c(cell_names, sampling_names)
				print(sprintf("%s's number is %d", cell, length(cell_names)))	
				#type_bycell[[cell]] <- subset(type_bycell[[cell]], cells = cell_names)
			}
		}
	}
	obj <- subset(obj, cells = sampling_names)
	return(obj)
}
