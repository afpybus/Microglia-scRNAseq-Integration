# PRELIMINARIES ################################

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#source("C:/Users/apybus3/Box/Wood Lab R Functions/WoodLabFunctions.R")
if (!require("pacman")) install.packages("pacman")
if (!require("BiocManager")) install.packages("BiocManager")
pacman::p_load(readxl,tidyverse,rio,Seurat,data.table,R.utils) 
pacman::p_load_gh("pengminshi/cFIT")

# 1 - Prep Arneson / TBI data set
# 2 - Analyze Arneson Monocytes
# 3 - Prep Sousa / LPS data set
# 4 - Integrate Data Sets
# 5 - cFIT




# 1 - PREP ARNESON DATA ##########

# Import the data
# data set is not included in GitHub repo due to file size limitations; download here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?&acc=GSE101901
# file GSE101901_DropSeqTBI.digital_expression.txt.gz, must decompress to .txt
df_TBI_in <- import("GEO Data Sets/DropSeqTBI.digital_expression.txt") 
genes_TBI = str_to_upper(df_TBI_in$V1)
dir.create("Figures/1 - Prep Arneson Data",recursive = TRUE,showWarnings = FALSE)

# Create count matrix
TBI_counts = df_TBI_in[,2:ncol(df_TBI_in)] %>% as.matrix()
rownames(TBI_counts)=genes_TBI

# Create data frame of metadata
meta_TBI = data.frame(Condition = case_when(
  str_detect(colnames(TBI_counts),"TBI") ~ "TBI",
  str_detect(colnames(TBI_counts),"Sham") ~ "Sham"
),
MouseID = case_when(
  str_detect(colnames(TBI_counts),"TBI1") ~ "TBI1",
  str_detect(colnames(TBI_counts),"TBI2") ~ "TBI2",
  str_detect(colnames(TBI_counts),"TBI3") ~ "TBI3",
  str_detect(colnames(TBI_counts),"Sham1") ~ "Sham1",
  str_detect(colnames(TBI_counts),"Sham2") ~ "Sham2",
  str_detect(colnames(TBI_counts),"Sham3") ~ "Sham3"
))
rownames(meta_TBI) = colnames(TBI_counts)

# Create Seurat object for use in Seurat package codes
TBI_SO <- CreateSeuratObject(counts = TBI_counts, project = "TBI",meta.data = meta_TBI) #, min.cells = 3, min.features = 200)

# Quality control: see distribution of number of features, number of counts  in all samples
TBI_SO[["percent.mt"]] <- PercentageFeatureSet(TBI_SO, pattern = "^MT-")
png("Figures/1 - Prep Arneson Data/MT RNA Violins.png",res=600,units="in",height=5,width=7)
VlnPlot(TBI_SO, features ="percent.mt")
dev.off()
png("Figures/1 - Prep Arneson Data/nFeature Violins.png",res=600,units="in",height=5,width=7)
VlnPlot(TBI_SO, features ="nFeature_RNA")
dev.off()
png("Figures/1 - Prep Arneson Data/nCount Violins.png",res=600,units="in",height=5,width=7)
VlnPlot(TBI_SO, features = "nCount_RNA")
dev.off()
png("Figures/1 - Prep Arneson Data/Features vs Counts Pre-Filter.png",res=600,units="in",height=5,width=6)
FeatureScatter(TBI_SO, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

# Filter out cells with too few (possibly droplets) or too many counts (possibly multiple cells) and cells with too much mito RNA
TBI_SO <- subset(TBI_SO, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 75)

png("Figures/1 - Prep Arneson Data/Features vs Counts Post-Filter.png",res=600,units="in",height=5,width=6)
FeatureScatter(TBI_SO, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

png("Figures/1 - Prep Arneson Data/MT RNA Violins Post-Filter.png",res=600,units="in",height=5,width=7)
VlnPlot(TBI_SO, features ="percent.mt")
dev.off()

# Normalize data with log transform
TBI_SO <- NormalizeData(TBI_SO)
TBI_SO <- FindVariableFeatures(TBI_SO, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(TBI_SO), 10)

# plot variable features 
png("Figures/1 - Prep Arneson Data/Variable Features.png",res=600,units="in",height=4,width=6)
VariableFeaturePlot(TBI_SO)
LabelPoints(plot = VariableFeaturePlot(TBI_SO), points = top10, repel = TRUE)
dev.off()

#Scale the data for PCA
all.genes <- rownames(TBI_SO)
TBI_SO <- ScaleData(TBI_SO, features = all.genes)

# Conduct PCA
TBI_SO <- RunPCA(TBI_SO, features = VariableFeatures(object = TBI_SO))

png("Figures/1 - Prep Arneson Data/PCA Loadings Charts PC1-2.png",res=600,units="in",height=4.5,width=6)
VizDimLoadings(TBI_SO, dims = 1:2, reduction = "pca")
dev.off()

png("Figures/1 - Prep Arneson Data/PCA Scores Plot PC1-2.png",res=600,units="in",height=4.5,width=6)
DimPlot(TBI_SO, reduction = "pca",group.by = "MouseID")
dev.off()

png("Figures/1 - Prep Arneson Data/PCA Heatmaps PC1-12.png",res=600,units="in",height=11.5,width=7)
DimHeatmap(TBI_SO, dims = 1:12, cells = 500, balanced = TRUE)
dev.off()

# # Determine appropriate dimension of the data by JackStraw analysis of PCs
# TBI_SO <- JackStraw(TBI_SO, num.replicate = 100, dims=30)
# TBI_SO <- ScoreJackStraw(TBI_SO, dims = 1:30)
# 
# png("Figures/1 - Prep Arneson Data/Jack Straw Plot.png",res=600,units="in",height=4,width=6)
# JackStrawPlot(TBI_SO, dims = 1:30)
# dev.off()

# PC Percent Variance Analysis
pct <- TBI_SO[["pca"]]@stdev / sum(TBI_SO[["pca"]]@stdev) * 100 # Determine percent of variation associated with each PC
cumu <- cumsum(pct)  # Calculate cumulative percents for each PC

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

png("Figures/1 - Prep Arneson Data/PCA Cumulative Variance.png",res=600,units="in",height=4,width=6)
plot(x=1:50,y=cumu,xlab="PC",ylab="Cumulative Variance (%)",pch=20)+
  abline(h=90,lty=2) +
  abline(v=co1,lty=2,col="red")
dev.off()

png("Figures/1 - Prep Arneson Data/PCA Elbow Plot.png",res=600,units="in",height=4,width=6)
ElbowPlot(TBI_SO,ndims = 50)
dev.off()

# Cluster the data using KNN of PCs and Louvain modularity optimization
TBI_SO <- FindNeighbors(TBI_SO, dims = 1:co1) %>%
  FindClusters(resolution = 0.5)

# Data visualization with nonlinear dim reduction UMAP
TBI_SO <- RunUMAP(TBI_SO, dims = 1:co1)

png("Figures/1 - Prep Arneson Data/UMAP Clusters.png",res=600,units="in",height=5,width=6)
DimPlot(TBI_SO, reduction = "umap",label=TRUE)
dev.off()

png("Figures/1 - Prep Arneson Data/UMAP by MouseID.png",res=600,units="in",height=5,width=6)
DimPlot(TBI_SO, reduction = "umap",group.by = "MouseID")
dev.off()

# Identify monocyte clusters

png("Figures/1 - Prep Arneson Data/Monocyte Genes by Cluster CTSS CD14.png",res=600,units="in",height=3.5,width=7)
VlnPlot(TBI_SO, features = c("CTSS", "CD14"))
dev.off()

png("Figures/1 - Prep Arneson Data/Monocyte Genes by Cluster AIF1 FCRLS.png",res=600,units="in",height=3.5,width=7)
VlnPlot(TBI_SO, features = c("AIF1","FCRLS"))
dev.off()

png("Figures/1 - Prep Arneson Data/Monocyte Genes by Cluster TMEM119 CCR5.png",res=600,units="in",height=3.5,width=7)
VlnPlot(TBI_SO, features = c("TMEM119","CCR5"))
dev.off()

png("Figures/1 - Prep Arneson Data/UMAP Monocyte Genes by Cluster.png",res=600,units="in",height=8,width=6)
FeaturePlot(TBI_SO, features = c("CTSS","CD14", "AIF1","FCRLS","TMEM119","CCR5"))
dev.off()


# Monocyte clusters are 4 and 14



# 2 - ANALYZE ARNESON MONOCYTES #############

dir.create("Figures/2 - Analyze Arneson Monocytes",recursive = TRUE,showWarnings = FALSE)

# Filter out all cells not in clusters 4 or 14
TBI_M <- TBI_SO[,which(TBI_SO$seurat_clusters %in% c(4,14))]

# Verify that count vs feature is linear
png("Figures/2 - Analyze Arneson Monocytes/Features vs Counts.png",res=600,units="in",height=5,width=6)
FeatureScatter(TBI_M, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

# Normalize data with log transform
TBI_M <- NormalizeData(TBI_M)
TBI_M <- FindVariableFeatures(TBI_M, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(TBI_M), 10)

# Plot variable features 
png("Figures/2 - Analyze Arneson Monocytes/Variable Features.png",res=600,units="in",height=4,width=6)
VariableFeaturePlot(TBI_M)
LabelPoints(plot = VariableFeaturePlot(TBI_M), points = top10, repel = TRUE)
dev.off()

# Scale the data for PCA
all.genes <- rownames(TBI_M)
TBI_M <- ScaleData(TBI_M, features = all.genes)

# Conduct PCA
TBI_M <- RunPCA(TBI_M, features = VariableFeatures(object = TBI_M))

png("Figures/2 - Analyze Arneson Monocytes/PCA Loadings Charts PC1-2.png",res=600,units="in",height=4.5,width=6)
VizDimLoadings(TBI_M, dims = 1:2, reduction = "pca")
dev.off()

png("Figures/2 - Analyze Arneson Monocytes/PCA Scores Plot PC1-2.png",res=600,units="in",height=4.5,width=6)
DimPlot(TBI_M, reduction = "pca",group.by = "MouseID")
dev.off()

png("Figures/2 - Analyze Arneson Monocytes/PCA Heatmaps PC1-12.png",res=600,units="in",height=11.5,width=7)
DimHeatmap(TBI_M, dims = 1:12, cells = 500, balanced = TRUE)
dev.off()

# # Determine appropriate dimension of the data by JackStraw analysis of PCs
# TBI_M <- JackStraw(TBI_M, num.replicate = 100, dims=30)
# TBI_M <- ScoreJackStraw(TBI_M, dims = 1:30)
# 
# png("Figures/2 - Analyze Arneson Monocytes/Jack Straw Plot.png",res=600,units="in",height=4,width=6)
# JackStrawPlot(TBI_M, dims = 1:30)
# dev.off()

# PC Percent Variance Analysis
pct <- TBI_M[["pca"]]@stdev / sum(TBI_M[["pca"]]@stdev) * 100 # Determine percent of variation associated with each PC
cumu <- cumsum(pct)  # Calculate cumulative percents for each PC

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

png("Figures/2 - Analyze Arneson Monocytes/PCA Cumulative Variance.png",res=600,units="in",height=4,width=6)
plot(x=1:50,y=cumu,xlab="PC",ylab="Cumulative Variance (%)",pch=20)+
  abline(h=90,lty=2) +
  abline(v=co1,lty=2,col="red")
dev.off()

png("Figures/2 - Analyze Arneson Monocytes/PCA Elbow Plot.png",res=600,units="in",height=4,width=6)
ElbowPlot(TBI_M,ndims = 50)
dev.off()

# Cluster the data using KNN of PCs and Louvain modularity optimization
TBI_M <- FindNeighbors(TBI_M, dims = 1:co1) %>%
  FindClusters(resolution = 0.5)

# Data visualization with nonlinear dim reduction UMAP
TBI_M <- RunUMAP(TBI_M, dims = 1:co1)

png("Figures/2 - Analyze Arneson Monocytes/UMAP Clusters.png",res=600,units="in",height=5,width=6)
DimPlot(TBI_M, reduction = "umap",label=TRUE)
dev.off()

# Find markers for the clusters
TBI_M.markers <- FindAllMarkers(TBI_M, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
TBI_M.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

png("Figures/2 - Analyze Arneson Monocytes/Monocyte Genes by Cluster 0.png",res=600,units="in",height=3.5,width=7)
VlnPlot(TBI_M, features = TBI_M.markers$gene[which(TBI_M.markers$cluster==0)[1:3]])
dev.off()

png("Figures/2 - Analyze Arneson Monocytes/Monocyte Genes by Cluster 1.png",res=600,units="in",height=3.5,width=7)
VlnPlot(TBI_M, features = TBI_M.markers$gene[which(TBI_M.markers$cluster==1)[1:3]])
dev.off()

png("Figures/2 - Analyze Arneson Monocytes/Monocyte Genes by Cluster 2.png",res=600,units="in",height=3.5,width=7)
VlnPlot(TBI_M, features = TBI_M.markers$gene[which(TBI_M.markers$cluster==2)[1:3]])
dev.off()

png("Figures/2 - Analyze Arneson Monocytes/Monocyte Genes by Cluster Custom.png",res=600,units="in",height=3.5,width=7)
VlnPlot(TBI_M, features = c("TMEM119","CXCL2","S100A6"))
dev.off()

png("Figures/2 - Analyze Arneson Monocytes/UMAP Monocyte Genes by Cluster.png",res=600,units="in",height=6,width=6)
FeaturePlot(TBI_M, features = c("TMEM119","CXCL2","S100A6","CD14"))
dev.off()

# Assign labels to clusters
new.cluster.ids <- c("Microglia", "High MT Content Cells", "Macrophage")
names(new.cluster.ids) <- levels(TBI_M)
TBI_M <- RenameIdents(TBI_M, new.cluster.ids)

png("Figures/2 - Analyze Arneson Monocytes/UMAP Clusters Labeled.png",res=600,units="in",height=5,width=6)
DimPlot(TBI_M, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()



# 3 - PREP SOUSA DATA ##############

dir.create("Figures/3 - Prep Sousa Data",recursive = TRUE,showWarnings = FALSE)

# Import the data
LPS <- fread("GEO Data Sets/GSM3182556_LPS.txt.gz")
CTRL <- fread("GEO Data Sets/GSM3182555_CTRL.txt.gz")

# Create count matrix
LPS_all = inner_join(LPS,CTRL,by="GENE")
LPS_counts = LPS_all[,2:ncol(LPS_all)] %>% as.matrix()
rownames(LPS_counts) = str_to_upper(LPS_all$GENE) 
colnames(LPS_counts) = c(
  str_c(rep("LPS",ncol(LPS)-1),c(1:{ncol(LPS)-1})),
  str_c(rep("CTRL",ncol(CTRL)-1),c(1:{ncol(CTRL)-1})))

# Create data frame of metadata
meta_LPS <- data.frame(Condition = case_when(
  str_detect(colnames(LPS_counts),"LPS") ~ "LPS",
  str_detect(colnames(LPS_counts),"CTRL") ~ "CTRL"
))
rownames(meta_LPS) = colnames(LPS_counts)



# Create Seurat object for use in Seurat package codes
LPS_SO <- CreateSeuratObject(counts = LPS_counts, project = "LPS",meta.data = meta_LPS) #, min.cells = 3, min.features = 200)

# Quality control: see distribution of number of features, number of counts  in all samples
png("Figures/3 - Prep Sousa Data/nFeature Violins.png",res=600,units="in",height=5,width=7)
VlnPlot(LPS_SO, features ="nFeature_RNA",group.by = "Condition")
dev.off()
png("Figures/3 - Prep Sousa Data/nCount Violins.png",res=600,units="in",height=5,width=7)
VlnPlot(LPS_SO, features = "nCount_RNA",group.by = "Condition")
dev.off()
png("Figures/3 - Prep Sousa Data/Features vs Counts Pre-Filter.png",res=600,units="in",height=5,width=6)
FeatureScatter(LPS_SO, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "Condition")
dev.off()

# Filter out cells with too few (possibly droplets) or too many counts (possibly multiple cells) 
LPS_SO <- subset(LPS_SO, subset = nFeature_RNA > 200 & nFeature_RNA < 4000)

png("Figures/3 - Prep Sousa Data/Features vs Counts Post-Filter.png",res=600,units="in",height=5,width=6)
FeatureScatter(LPS_SO, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by="Condition")
dev.off()

# Perform Integration of LPS and CTRL data
# split the dataset into a list of two seurat objects 
LPS.list <- SplitObject(LPS_SO, split.by = "Condition")

# normalize and identify variable features for each dataset independently
LPS.list <- lapply(X = LPS.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = LPS.list)

anchors <- FindIntegrationAnchors(object.list = LPS.list, anchor.features = features)
# this command creates an 'integrated' data assay
LPS.combined <- IntegrateData(anchorset = anchors)
DefaultAssay(LPS.combined) <- "integrated"

# Find variable features
LPS.combined <- FindVariableFeatures(LPS.combined, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(LPS.combined), 10)

# plot variable features 
png("Figures/3 - Prep Sousa Data/Variable Features.png",res=600,units="in",height=4,width=6)
VariableFeaturePlot(LPS.combined)
LabelPoints(plot = VariableFeaturePlot(LPS.combined), points = top10, repel = TRUE)
dev.off()

#Scale the data for PCA
all.genes <- rownames(LPS.combined)
LPS.combined <- ScaleData(LPS.combined, features = all.genes)

# Conduct PCA
LPS.combined <- RunPCA(LPS.combined, features = VariableFeatures(object = LPS.combined))

png("Figures/3 - Prep Sousa Data/PCA Loadings Charts PC1-2.png",res=600,units="in",height=4.5,width=6)
VizDimLoadings(LPS.combined, dims = 1:2, reduction = "pca")
dev.off()

png("Figures/3 - Prep Sousa Data/PCA Scores Plot PC1-2.png",res=600,units="in",height=4.5,width=6)
DimPlot(LPS.combined, reduction = "pca",group.by = "Condition")
dev.off()

png("Figures/3 - Prep Sousa Data/PCA Heatmaps PC1-12.png",res=600,units="in",height=11.5,width=7)
DimHeatmap(LPS.combined, dims = 1:12, cells = 500, balanced = TRUE)
dev.off()

# # Determine appropriate dimension of the data by JackStraw analysis of PCs
# LPS.combined <- JackStraw(LPS.combined, num.replicate = 100, dims=30)
# LPS.combined <- ScoreJackStraw(LPS.combined, dims = 1:30)
# 
# png("Figures/3 - Prep Sousa Data/Jack Straw Plot.png",res=600,units="in",height=4,width=6)
# JackStrawPlot(LPS.combined, dims = 1:30)
# dev.off()

# PC Percent Variance Analysis
pct <- LPS.combined[["pca"]]@stdev / sum(LPS.combined[["pca"]]@stdev) * 100 # Determine percent of variation associated with each PC
cumu <- cumsum(pct)  # Calculate cumulative percents for each PC

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

png("Figures/3 - Prep Sousa Data/PCA Cumulative Variance.png",res=600,units="in",height=4,width=6)
plot(x=1:50,y=cumu,xlab="PC",ylab="Cumulative Variance (%)",pch=20)+
  abline(h=90,lty=2) +
  abline(v=co1,lty=2,col="red")
dev.off()

png("Figures/3 - Prep Sousa Data/PCA Elbow Plot.png",res=600,units="in",height=4,width=6)
ElbowPlot(LPS.combined,ndims = 50)
dev.off()

# Cluster the data using KNN of PCs and Louvain modularity optimization
LPS.combined <- FindNeighbors(LPS.combined, dims = 1:co1) %>%
  FindClusters(resolution = 0.5)

# Data visualization with nonlinear dim reduction UMAP
LPS.combined <- RunUMAP(LPS.combined, dims = 1:co1)

png("Figures/3 - Prep Sousa Data/UMAP Clusters.png",res=600,units="in",height=5,width=6)
DimPlot(LPS.combined, reduction = "umap",label=TRUE)
dev.off()

png("Figures/3 - Prep Sousa Data/UMAP by Condition.png",res=600,units="in",height=5,width=6)
DimPlot(LPS.combined, reduction = "umap",group.by = "Condition")
dev.off()

# Identifyclusters

# Find markers for the clusters
LPS.markers <- FindAllMarkers(LPS.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
LPS.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

for(i in 0:{length(unique(LPS.combined$seurat_clusters))-1}){
  png(paste0("Figures/3 - Prep Sousa Data/Marker Genes by Cluster ",i," .png"),res=600,units="in",height=3.5,width=7)
  print({VlnPlot(LPS.combined, features = LPS.markers$gene[which(LPS.markers$cluster==i)[1:3]])})
  dev.off()
}

png("Figures/3 - Prep Sousa Data/UMAP Monocyte Genes by Cluster.png",res=600,units="in",height=8,width=6)
FeaturePlot(LPS.combined, features = c("CTSS","CD14", "AIF1","FCRLS","TMEM119","CCR5"))
dev.off()



# 4 - INTEGRATE DATA SETS ###############

dir.create("Figures/4 - Integrate Data Sets",recursive = TRUE,showWarnings = FALSE)

# Perform Integration of LPS and CTRL data
combined_SO <- merge(TBI_M, y = LPS.combined, add.cell.ids = c("Arnseon", "Sousa"), project = "Integration")

# split the dataset into a list of two seurat objects 
SO.list <- SplitObject(combined_SO, split.by = "Condition")

# normalize and identify variable features for each dataset independently
SO.list <- lapply(X = SO.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = SO.list)
anchors <- FindIntegrationAnchors(object.list = SO.list, anchor.features = features)

# this command creates an 'integrated' data assay
SO.combined <- IntegrateData(anchorset = anchors)
DefaultAssay(SO.combined) <- "integrated"

# Find variable features
SO.combined <- FindVariableFeatures(SO.combined, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(SO.combined), 10)

# plot variable features 
png("Figures/4 - Integrate Data Sets/Variable Features.png",res=600,units="in",height=4,width=6)
VariableFeaturePlot(SO.combined)
LabelPoints(plot = VariableFeaturePlot(SO.combined), points = top10, repel = TRUE)
dev.off()

#Scale the data for PCA
all.genes <- rownames(SO.combined)
SO.combined <- ScaleData(SO.combined, features = all.genes)

# Conduct PCA
SO.combined <- RunPCA(SO.combined, features = VariableFeatures(object = SO.combined))

png("Figures/4 - Integrate Data Sets/PCA Loadings Charts PC1-2.png",res=600,units="in",height=4.5,width=6)
VizDimLoadings(SO.combined, dims = 1:2, reduction = "pca")
dev.off()

png("Figures/4 - Integrate Data Sets/PCA Scores Plot PC1-2.png",res=600,units="in",height=4.5,width=6)
DimPlot(SO.combined, reduction = "pca",group.by = "Condition")
dev.off()

png("Figures/4 - Integrate Data Sets/PCA Heatmaps PC1-12.png",res=600,units="in",height=11.5,width=7)
DimHeatmap(SO.combined, dims = 1:12, cells = 500, balanced = TRUE)
dev.off()

# # Determine appropriate dimension of the data by JackStraw analysis of PCs
# SO.combined <- JackStraw(SO.combined, num.replicate = 100, dims=30)
# SO.combined <- ScoreJackStraw(SO.combined, dims = 1:30)
# 
# png("Figures/4 - Integrate Data Sets/Jack Straw Plot.png",res=600,units="in",height=4,width=6)
# JackStrawPlot(SO.combined, dims = 1:30)
# dev.off()

# PC Percent Variance Analysis
pct <- SO.combined[["pca"]]@stdev / sum(SO.combined[["pca"]]@stdev) * 100 # Determine percent of variation associated with each PC
cumu <- cumsum(pct)  # Calculate cumulative percents for each PC

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

png("Figures/4 - Integrate Data Sets/PCA Cumulative Variance.png",res=600,units="in",height=4,width=6)
plot(x=1:50,y=cumu,xlab="PC",ylab="Cumulative Variance (%)",pch=20)+
  abline(h=90,lty=2) +
  abline(v=co1,lty=2,col="red")
dev.off()

png("Figures/4 - Integrate Data Sets/PCA Elbow Plot.png",res=600,units="in",height=4,width=6)
ElbowPlot(SO.combined,ndims = 50)
dev.off()

# Cluster the data using KNN of PCs and Louvain modularity optimization
SO.combined <- FindNeighbors(SO.combined, dims = 1:co1) %>%
  FindClusters(resolution = 0.5)

# Data visualization with nonlinear dim reduction UMAP
SO.combined <- RunUMAP(SO.combined, dims = 1:co1)

png("Figures/4 - Integrate Data Sets/UMAP Clusters.png",res=600,units="in",height=5,width=6)
DimPlot(SO.combined, reduction = "umap",label=TRUE)
dev.off()

png("Figures/4 - Integrate Data Sets/UMAP by Condition.png",res=600,units="in",height=5,width=6)
DimPlot(SO.combined, reduction = "umap",group.by = "Condition")
dev.off()

# Identify clusters
# For performing differential expression after integration, we switch back to the original data
DefaultAssay(SO.combined) <- "RNA"

# Find markers for the clusters
SO.markers <- FindAllMarkers(SO.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
SO.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

for(i in 0:{length(unique(SO.combined$seurat_clusters))-1}){
  png(paste0("Figures/4 - Integrate Data Sets/Marker Genes by Cluster ",i," .png"),res=600,units="in",height=3.5,width=7)
  print({VlnPlot(SO.combined, features = SO.markers$gene[which(SO.markers$cluster==i)[1:3]])})
  dev.off()
  
  gene_export = as.matrix(SO.markers$gene[which(SO.markers$cluster==i)])
  export(gene_export,file=paste0("Figures/4 - Integrate Data Sets/Top Genes Cluster ",i,".txt"),)
}

# png("Figures/4 - Integrate Data Sets/UMAP Monocyte Genes by Cluster.png",res=600,units="in",height=8,width=6)
# FeaturePlot(SO.combined, features = c("CTSS","CD14", "AIF1","FCRLS","TMEM119","CCR5"))
# dev.off()

png(paste0("Figures/4 - Integrate Data Sets/MG-Specific Genes by Cluster.png"),res=600,units="in",height=3,width=7)
VlnPlot(SO.combined, features = c("TMEM119","CCR5"))
dev.off()

VlnPlot(SO.combined, features = c("TMEM119","CCL5"))


genes = rownames(combined_SO@assays$RNA)
"H2EA" %in% genes
genes[which(str_detect(genes,"H2"))]



# 5 - cFIT ############

dir.create("Figures/5 - cFIT",recursive = TRUE,showWarnings = FALSE)

data.list = split_dataset_by_batch(X=t(as.matrix(SO.combined@assays$RNA@counts)), 
                                   batch = SO.combined@meta.data$Condition, 
                                   labels = SO.combined@meta.data$seurat_clusters, 
                                   metadata = SO.combined@meta.data, 
                                   dataset.name = 'Integrated')

# select 2000 highly variable genes
genes = select_genes(data.list$X.list, ngenes=2000, verbose=F)

# data preprocessing
exprs.list = preprocess_for_integration(data.list$X.list, genes, scale.factor=10^4, scale=F, center=F)

# run cFIT
int.out = CFITIntegrate(X.list=exprs.list, r=15, verbose=F, max.niter = 100, seed=0)

# save the results
saveRDS(int.out, file='Figures/5 - cFIT/cFIT Arneson Sousa.rds')

# Obtain the ncell-by-ngene expression matrix
exprs.int = do.call(rbind, int.out$H.list) %*% t(int.out$W)

# Obstain the ncell-by-r low dimensional representation
Hnorm = do.call(rbind, int.out$H.list) %*% diag(colSums(int.out$W))


# Visualize the data with UMAP

umap.out = plot_umap(X=Hnorm,  
                     pca = NULL, n_components = 2, n_neighbors = 50, min_dist = 0.1, # umap parameters
                     point.size = 0.6, alpha=0.8, title=NULL, legend.name='cell type', # figure parameters
                     seed=42)

p1 = plot_umap(labels=SO.combined@meta.data$Condition, point.size = 0.5, alpha=0.5, legend.name='Condition', emb=umap.out$emb)$p 
p2 = plot_umap(labels=SO.combined@meta.data$seurat_clusters, point.size = 0.5, alpha=0.5, legend.name='Cluster', emb=umap.out$emb)$p

png("Figures/5 - cFIT/cFIT Integration by Condition.png",res=600,units="in",height=3,width=4)
p1
dev.off()

png("Figures/5 - cFIT/cFIT Integration by Cluster.png",res=600,units="in",height=3,width=4)
p2
dev.off()



####################








