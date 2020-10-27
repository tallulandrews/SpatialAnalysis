require(Seurat)
require(ggplot2)


slice <- "PSC011_4_D1"
set.seed(101)

obj <- Load10X_Spatial(paste("./", slice, sep=""), 
		filename="filtered_feature_bc_matrix.h5", 
		slice=slice)

obj@meta.data$orig.ident <- rep(slice, ncol(obj))

# QC
plot_count <- SpatialFeaturePlot(obj, features = "nCount_Spatial") + 
		theme(legend.position = "right")
plot_feature <- SpatialFeaturePlot(obj, features = "nFeature_Spatial") + 
		theme(legend.position = "right")
plot_count + plot_feature


obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

plot_mt <- SpatialFeaturePlot(obj, features = "percent.mt") + 
		theme(legend.position = "right")
plot_mt


# FeatureSelection
# Filter
detect_freq <- Matrix::rowMeans(obj@assays$Spatial@counts > 0)
highest <- apply(obj@assays$Spatial@counts, 1, max)

obj <- obj[detect_freq > 0.05 & highest >= 3,]
obj <- obj[!grepl("^MT-", rownames(obj)),]


# Normalization
obj <- SCTransform(obj, assay = "Spatial", verbose = FALSE)

mean_expr <- Matrix::rowMeans(obj@assays$Spatial@data)

# Pipeline
obj <- RunPCA(obj, assay = "SCT", verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:20)
obj <- FindClusters(obj, verbose = FALSE)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:20)

p1 <- DimPlot(obj, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(obj, label = TRUE, label.size = 3)

p2

# spatial plot of pca scores

obj@meta.data$PC1_score <- obj@reductions$pca@cell.embeddings[,1]
obj@meta.data$PC2_score <- obj@reductions$pca@cell.embeddings[,2]
obj@meta.data$PC3_score <- obj@reductions$pca@cell.embeddings[,3]
obj@meta.data$PC4_score <- obj@reductions$pca@cell.embeddings[,4]
obj@meta.data$PC5_score <- obj@reductions$pca@cell.embeddings[,5]
obj@meta.data$PC6_score <- obj@reductions$pca@cell.embeddings[,6]
pca_plot1 <- SpatialFeaturePlot(obj, features = c("PC1_score","PC2_score", "PC3_score", "PC4_score", "PC5_score", "PC6_score"))


obj@meta.data$PC7_score <- obj@reductions$pca@cell.embeddings[,7]
obj@meta.data$PC8_score <- obj@reductions$pca@cell.embeddings[,8]
obj@meta.data$PC9_score <- obj@reductions$pca@cell.embeddings[,9]
obj@meta.data$PC10_score <- obj@reductions$pca@cell.embeddings[,10]
obj@meta.data$PC11_score <- obj@reductions$pca@cell.embeddings[,11]
obj@meta.data$PC12_score <- obj@reductions$pca@cell.embeddings[,12]
pca_plot2 <- SpatialFeaturePlot(obj, features = c("PC7_score","PC8_score", "PC9_score", "PC10_score", "PC11_score", "PC12_score"))


require(fgsea)

pathway_files <- c("C:/Users/tandrews/Documents/UHNSonya/ExternalData/MSigDbHalmarkPathways.gmt",
			"C:/Users/tandrews/Documents/UHNSonya/ExternalData/MSigDb_immune_signatures_c7.all.v7.1.symbols.gmt",
			"C:/Users/tandrews/Documents/UHNSonya/ExternalData/MSigDb_curated_c2.all.v7.1.symbols.gmt"
			)

immune_pathways <- fgsea::gmtPathways(pathway_files[2])
hallmark_pathways <- fgsea::gmtPathways(pathway_files[1])
curated_pathways <- fgsea::gmtPathways(pathway_files[3])

source("C:/Users/tandrews/Documents/UHNSonya/AutoAnnotation/Setup_autoannotation.R")

PCA_Top_Pos <- c();
PCA_Top_Neg <- c();

auto_anno_marks_as_pathways <- list()
types <- unique(map1_markers$label)
types <- types[types != "None"]
for (type in types){
	genes <- as.character(map1_markers[map1_markers[,2] == type,1])
	genes <- genes[!grepl("^RPS", genes)]
	genes <- genes[!grepl("^RPL", genes)]
	auto_anno_marks_as_pathways[[type]] <- genes
}
for (i in 1:ncol(obj@reductions$pca@feature.loadings)) {
	print(i)
	out_types <- fgsea(auto_anno_marks_as_pathways, obj@reductions$pca@feature.loadings[,i])
	print(out_types[out_types$padj == min(out_types$padj),])

	tmp <- sort(obj@reductions$pca@feature.loadings[,i])
	PCA_Top_Pos <- cbind(PCA_Top_Pos, tail(names(tmp), 25))
	PCA_Top_Neg <- cbind(PCA_Top_Neg, head(names(tmp), 25))
}

colnames(PCA_Top_Pos) <- paste("PC", 1:ncol(PCA_Top_Pos))
write.table(PCA_Top_Pos, file=paste(slice, "Top_Pos_PCA_genes.txt",sep="_"), sep="\t", quote=TRUE, row.names=F, col.names=T)

colnames(PCA_Top_Neg) <- paste("PC", 1:ncol(PCA_Top_Neg))
write.table(PCA_Top_Neg, file=paste(slice, "Top_Neg_PCA_genes.txt",sep="_"), sep="\t", quote=TRUE, row.names=F, col.names=T)



#out_hallmark <- fgsea(hallmark_pathways, obj@reductions$pca@feature.loadings[,1], nperm=1000000)
#out_curated <- fgsea(curated_pathways, obj@reductions$pca@feature.loadings[,1], nperm=1000000)
#out_immune <- fgsea(immune_pathways, obj@reductions$pca@feature.loadings[,1], nperm=1000000)

#component <- 8; obj@meta.data$PC_score <- obj@reductions$pca@cell.embeddings[,component]; SpatialFeaturePlot(obj, features = "PC_score")

#### specific key genes

marks <- SpatialFeaturePlot(obj, features = c("HBA1", "CYP3A4", "ALB", "IGKC", "MARCO", "GLUL"))

# HVGs
obj <- FindVariableFeatures(obj, assay="SCT", nfeatures=5000)
obj <- ScaleData(obj)

obj <- FindSpatiallyVariableFeatures(obj, assay="SCT", features = VariableFeatures(obj)[1:5000], 
    selection.method = "markvariogram")
spatial_features <- SpatiallyVariableFeatures(obj, selection.method = "markvariogram")

sp1 <- SpatialFeaturePlot(obj, features = spatial_features[1:25])
sp2 <- SpatialFeaturePlot(obj, features = VariableFeatures(obj)[1:25])

saveRDS(obj, paste(slice, "Processed_Spatial_Seurat.rds", sep="_"))











###### Cell-type Markers ########

liver_m <- read.table("Liver_markers.csv", sep=",", header=T)
immune_m <- read.table("Immune_markers.csv", sep=",", header=F)
general_m <- read.table("SoupXGenesets_markers.csv", sep=",", header=T)

tmp <- mean_expr[immune_m[,1]]
tmp <- names(tmp)[tmp>0.5 & !is.na(tmp)]
immune_m[immune_m[,1] %in% tmp,]

tmp <- mean_expr[liver_m[,1]]
tmp <- names(tmp)[tmp>0.5 & !is.na(tmp)]
liver_m[liver_m[,1] %in% tmp,]

SpatialFeaturePlot(obj, features = c("F10", "ENG", "FCN3", "GLUL", "ABCC3", "KRT18"))

SpatialFeaturePlot(obj, features = c("MARCO", "IGKC", "CD163", "HLA-DPA1"))



SpatialFeaturePlot(obj, features = 
	c("C1QC", "CD74", "CTSB", "FTH1", "FTL", "IGLC2", "IGKC", 
		"PDK4", "S100A8", "S100A9", "SDC1", "SLC40A1", "SSR4", "TUBB",
		"TUBA1B"))


good_markers <- c("MARCO", "CYP3A4", "CYP2E1", "IGKC", "GLUL", "CD163", "ALB", "CDC3D")


c("ALB, "CYP2E1", "GLUL", "IGKC", "CD163")

marks_mac <- c("C1QC", "AIF1l", "MARCO", "BACH1", "CCL2", "CCR2")

SpatialFeaturePlot(obj, features = good_markers)


SpatialFeaturePlot(obj, features = 
	c("ALB", "PTPRC", "MARCO", "CYP2E1", "GLUL", "IGKC"))