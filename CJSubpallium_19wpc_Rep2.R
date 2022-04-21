library(Seurat)
library(dplyr)
library(ggplot2)

#Seurat file for 19wk Marmoset fetal brain Subpallium replicate 2
# Make Seurat object using filtered feature matrices.

CJSubpallium19wkRep2_dir <- '~/MarmosetWork2022/CountJobs/CJSubpallium19wkRep2/outs/filtered_feature_bc_matrix'
CJSubpallium19wkRep2.data <- Read10X(data.dir = CJSubpallium19wkRep2_dir)
CJSubpallium19wkRep2 <- CreateSeuratObject(counts = CJSubpallium19wkRep2.data, project = "CJSubpallium19wkRep2", min.cells = 3, min.features = 200)
CJSubpallium19wkRep2[["percent.mt"]] <- PercentageFeatureSet(CJSubpallium19wkRep2, pattern = "^MT-")
Vlnplot1 <- VlnPlot(CJSubpallium19wkRep2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("~/MarmosetWork2022/CJSubpallium19wkRep2DataAnalysis/VlnplotQC1.jpeg", width = 12, height = 7, units = c("in"), dpi = 300)
dev.off()

plot1 <- FeatureScatter(CJSubpallium19wkRep2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(CJSubpallium19wkRep2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plotcombined <- plot1 + plot2
ggsave("~/MarmosetWork2022/CJSubpallium19wkRep2DataAnalysis/featureRNAplot.jpeg", width = 16, height = 5, units = c("in"), dpi = 300)
dev.off()

#based on the QC, choosing cells above 500 and less that 8000
CJSubpallium19wkRep2 <- subset(CJSubpallium19wkRep2, subset = nFeature_RNA > 500 & nFeature_RNA < 8000)

#Normalize
CJSubpallium19wkRep2 <- NormalizeData(CJSubpallium19wkRep2, normalization.method = "LogNormalize", scale.factor = 10000)

#Find highly variable features
CJSubpallium19wkRep2 <- FindVariableFeatures(CJSubpallium19wkRep2, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(CJSubpallium19wkRep2), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(CJSubpallium19wkRep2)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plotcombinedvariance <- plot1 + plot2
ggsave("~/MarmosetWork2022/CJSubpallium19wkRep2DataAnalysis/Varianceplot.jpeg", width = 14, height =6, units = c("in"), dpi = 300)
dev.off()

#Scaling the data 
all.genes <- rownames(CJSubpallium19wkRep2)
CJSubpallium19wkRep2<- ScaleData(CJSubpallium19wkRep2, features = all.genes)

#Perform linear dimensional reduction
CJSubpallium19wkRep2 <- RunPCA(CJSubpallium19wkRep2, features = VariableFeatures(object = CJSubpallium19wkRep2))

#determine the dimensionality
CJSubpallium19wkRep2 <- JackStraw(CJSubpallium19wkRep2, num.replicate = 100)
CJSubpallium19wkRep2<- ScoreJackStraw(CJSubpallium19wkRep2, dims = 1:20)
JackStrawPlot(CJSubpallium19wkRep2, dims = 1:15)

#Elbow plot for finding the dimensionality
ElbowPlot(CJSubpallium19wkRep2)

#Cluster the cells based on the PCs 
#In this case choosing 19 based on Elbow Plot
CJSubpallium19wkRep2 <- FindNeighbors(CJSubpallium19wkRep2, dims = 1:19)
CJSubpallium19wkRep2 <- FindClusters(CJSubpallium19wkRep2, resolution = 0.5)
CJSubpallium19wkRep2 <- RunUMAP(CJSubpallium19wkRep2, dims = 1:19)
DimPlot(CJSubpallium19wkRep2, reduction = "umap")
ggsave("~/MarmosetWork2022/CJSubpallium19wkRep2DataAnalysis/UMAP.jpeg", width = 7, height = 5, units = c("in"), dpi = 300)
dev.off()

#Save the Seurat file
#This can be easily loaded for cell type identification
saveRDS(CJSubpallium19wkRep2, file = "~/MarmosetWork2022/ScriptsMarmosetWork2022/CJSubpallium19wkRep2Seur.rds")

#Find the markers
markers <- FindAllMarkers(CJSubpallium19wkRep2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(markers, file ="~/MarmosetWork2022/CJSubpallium19wkRep2DataAnalysis/markersCJCortex19wkRep2.txt", sep = "\t", quote=F,row.names =T)
