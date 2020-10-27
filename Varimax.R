files <- c("C73_A1/C73_A1_Processed_Spatial_Seurat.rds",
		"C73_B1/C73_B1_Processed_Spatial_Seurat.rds",
		"C73_C1/C73_C1_Processed_Spatial_Seurat.rds",
		"C73_D1/C73_D1_Processed_Spatial_Seurat.rds",
		"PSC011_4_A1/PSC011_4_A1_Processed_Spatial_Seurat.rds",
		"PSC011_4_B1/PSC011_4_B1_Processed_Spatial_Seurat.rds",
		"PSC011_4_C1/PSC011_4_C1_Processed_Spatial_Seurat.rds",
		"PSC011_4_D1/PSC011_4_D1_Processed_Spatial_Seurat.rds")

source("RBO_from_gespeR.R")

#### From Delaram ####
get_varimax_rotated <- function(gene_exp_matrix, loading_matrix){
  
  ## gene_exp_matrix: gene expression matrix. rows named as genes and columns named as UMIs.
  ##                  attention: Data SHOULD be SCALED.
  ##                  Potentially gained from Seurat GetAssayData() function
  ## loading_matrix: the PCA loading matrix with rows named as gene names and 
  ##                 columns named as PC numbers. Potentially gained from Seurat Loadings() function
  
  ######### Varimax rotation
  initial_data <- t(gene_exp_matrix[rownames(gene_exp_matrix) %in% rownames(loading_matrix),])
  
  ## apply varimax rotation on the loadings
  varimax_res <- varimax(loading_matrix)
  rotatedLoadings <- varimax_res$loadings
  class(rotatedLoadings) <- "matrix";
  ## calculating the PC scores matrix
  invLoadings     <- t(pracma::pinv(rotatedLoadings))
  scores          <- scale(initial_data) %*% invLoadings ## this second scaling is not necessary
  
  ## compacting the rotated loading and score matrices in a list
  rotated_data <- list(rotLoadings=rotatedLoadings, rotScores = scores)
  return(rotated_data)
}



gene_lists <- c()
gene_lists_names <- c();

gene_lists_del <- c()
gene_lists_names_del <- c();


npcs <- 12;

gene_loadings_matrix <- c();

all_vmax_components <- list();

for (f in files) {
	tag <- unlist(strsplit(f, "\\/"))[1]
	print(tag)
	obj <- readRDS(f)
	outfile <- sub("Processed", "varimax", f)
	require(Seurat)

	# Delaram Varimax
#	set.seed(1029)
#	rotated <- get_varimax_rotated(gene_exp_matrix=as.matrix(obj@assays$Spatial@data), 
#						loading_matrix=obj@reductions$pca@feature.loadings[,1:npcs])
#	varimax_pca <- CreateDimReducObject(
#				embeddings = rotated$rotScores, 
#				loadings=rotated$rotLoadings, 
#				key = "VMD_", assay = DefaultAssay(obj))
#
#	obj@reductions[["varimax_del"]] <- varimax_pca
#	obj <- ProjectDim(obj, reduction="varimax_del")
#
#	colnames(rotated$rotLoadings) <- paste("RotPC_del", 1:ncol(rotated$rotLoadings), sep="_")
#	colnames(rotated$rotScores) <- paste("RotPC_del", 1:ncol(rotated$rotScores), sep="_")
#
#
#	obj@meta.data <- cbind(obj@meta.data, rotated$rotScores[,1:npcs])
#
#	png(paste(tag, "varimax_del_plot.png", sep="_"), width=10, height=7, unit="in", res=300)
#	print(SpatialFeaturePlot(obj, features=c("RotPC_del_1", "RotPC_del_2", "RotPC_del_3", "RotPC_del_4", "RotPC_del_5", "RotPC_del_6")))
#	dev.off()
#
#	png(paste(tag, "varimax_Del_plot2.png", sep="_"), width=10, height=7, unit="in", res=300)
#	print(SpatialFeaturePlot(obj, features=c("RotPC_del_7", "RotPC_del_8", "RotPC_del_9", "RotPC_del_10", "RotPC_del_11", "RotPC_del_12")))
#	dev.off()
#
#	for (i in 1:npcs) {
#		scores <- rotated$rotLoadings[,i]
#		scores <- sort(scores)
#		genes <- c(head(names(scores),50), tail(names(scores), 50))
#		gene_lists_del <- cbind(gene_lists_del, as.character(genes))
#		gene_lists_names_del <- c(gene_lists_names_del, paste(tag, i, sep="_"))
#	}



	# Old stuff
	set.seed(1029)
	pca_loading <- obj@reductions$pca@feature.loadings[,1:npcs]
	pca_loading <- pca_loading/rowSums(pca_loading)

	set.seed(1029)
	test_pca <- varimax(pca_loading)
	dim(test_pca$rotmat)

	require(scales)

	rotated_pca <- obj@reductions$pca@cell.embeddings[,1:npcs] %*% test_pca$rotmat
	loading_mat <- obj@reductions$pca@feature.loadings[,1:npcs] %*% test_pca$rotmat
	colnames(loading_mat) <- paste(tag, "RotPC", 1:ncol(loading_mat), sep="_")

	if (length(gene_loadings_matrix) == 0) {
		gene_loadings_matrix <- loading_mat
	} else {
		all_genes <- sort(unique(c(rownames(gene_loadings_matrix), rownames(loading_mat))))
		gene_loadings_matrix <- gene_loadings_matrix[match(all_genes, rownames(gene_loadings_matrix)),]
		loading_mat <- loading_mat[match(all_genes, rownames(loading_mat)),]
		loading_mat[is.na(loading_mat)] <- 0
		gene_loadings_matrix[is.na(gene_loadings_matrix)] <- 0;
		rownames(gene_loadings_matrix) <- all_genes
		rownames(loading_mat) <- all_genes
		gene_loadings_matrix <- cbind(gene_loadings_matrix, loading_mat)
	}


	colnames(rotated_pca) <- paste("RotPC", 1:ncol(rotated_pca), sep="_")

	obj@meta.data <- cbind(obj@meta.data, rotated_pca[,1:npcs])

	png(paste(tag, "varimax_plot.png", sep="_"), width=10, height=7, unit="in", res=300)
	print(SpatialFeaturePlot(obj, features=c("RotPC_1", "RotPC_2", "RotPC_3", "RotPC_4", "RotPC_5", "RotPC_6")))
	dev.off()

	png(paste(tag, "varimax_plot2.png", sep="_"), width=10, height=7, unit="in", res=300)
	print(SpatialFeaturePlot(obj, features=c("RotPC_7", "RotPC_8", "RotPC_9", "RotPC_10", "RotPC_11", "RotPC_12")))
	dev.off()



	varimax_pca <- CreateDimReducObject(
				embeddings = rotated_pca, 
				loadings=loading_mat, 
				key = "VM_", assay = DefaultAssay(obj))

	obj@reductions[["varimax"]] <- varimax_pca
	obj <- ProjectDim(obj, reduction="varimax")

	saveRDS(obj, outfile)

	for (i in 1:npcs) {
		scores <- loading_mat[,i]
		scores <- sort(scores)
		genes <- c(head(names(scores),50), tail(names(scores), 50))
		gene_lists <- cbind(gene_lists, as.character(genes))
		gene_lists_names <- c(gene_lists_names, paste(tag, i, sep="_"))
	}
	all_vmax_components[[tag]] <- obj@reductions$varimax@feature.loadings.projected
}


## Cluster components:
set.seed(1927)

gene_lists_old <- gene_lists
gene_lists_names_old <- gene_lists_names

gene_lists <- gene_lists_old
gene_lists_names <- gene_lists_names_old

require(proxy)

colnames(gene_lists) <- gene_lists_names

pairs <- t(combn(1:ncol(gene_lists), 2))
similarity <- apply(pairs, 1, function(x) {length(intersect(gene_lists[,x[1]], gene_lists[,x[2]]))})

sim_mat <- matrix(0, nrow=ncol(gene_lists), ncol=ncol(gene_lists))
for (i in 1:nrow(pairs)) {
	sim_mat[pairs[i,1], pairs[i,2]] <- similarity[i]
}

head(pairs[order(similarity, decreasing=T),], 20)

diff <- matrix(100, nrow=ncol(gene_lists), ncol=ncol(gene_lists))
diff <- diff-sim_mat
diff <- diff - t(sim_mat)
diag(diff) <- 0

require("gplots")
hmap <- heatmap.2(diff, scale="none", trace="none", col=colorRampPalette(c("black", "white"))(20))

clusters <- cutree(as.hclust(hmap$rowDendrogram), h = 110) 

sort(table(clusters))
colnames(gene_lists)[clusters==8]

genes <- c(unlist(gene_lists[1:50,clusters==46]))

genes1 <- c(unlist(gene_lists[1:50, colnames(gene_lists) %in% c("PSC011_4_A1_2", "PSC011_4_B1_2", "PSC011_4_C1_2", "PSC011_4_D3")]))#,
#		unlist(gene_lists[51:100, colnames(gene_lists) %in% c("PSC011_4_D1_6")]))

genes2 <- c(unlist(gene_lists[51:100, colnames(gene_lists) %in% c("PSC011_4_A1_4", "PSC011_4_B1_4", "PSC011_4_C1_4")]),
		unlist(gene_lists[1:50, colnames(gene_lists) %in% c("PSC011_4_D1_6")]))



