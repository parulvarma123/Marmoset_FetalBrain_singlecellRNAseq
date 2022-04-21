library(Seurat)
library(dplyr)
library(ggplot2)

#Load the previously saved Seurat file for CJ Cerebral Cortex
CJ19wkCortexRep2 <- readRDS(file = '~/MarmosetWork2022/ScriptsMarmosetWork2022/CJCortex19wkRep2Seur.rds')

#Identify the clusters based on the gene expression profiles using FindAllMarkers function of Seurat
#Renaming the clusters after marker analysis for each cluster

CJ19wkCortexRep2 <- RenameIdents(CJ19wkCortexRep2,
                                 '0' = "Radial Glial Cells", 
                                 '1' = "Oligodendrocyte Progenitor Cells", 
                                 '2' = "Astrocytes",
                                 '3' = "Unknown", 
                                 '4' = "Deep Layer Neurons", 
                                 '5' = "Unknown", 
                                 '6' = "Ventral Intermediate Progenitors", 
                                 '7' = "Radial Glial Cells", 
                                 '8' = "Microglia", 
                                 '9' = "Deep Layer Neurons",
                                 '10' = "Ventral Intermediate Progenitors", 
                                 '11' = "Pericytes", 
                                 '12' = "Developing Oligodendrocytes", 
                                 '13' = "Endothelial Cells", 
                                 '14' = "Dividing Progenitors",
                                 '15' = "Ventral Intermediate Progenitors",
                                 '16' = "Microglia",
                                 '17' = "Extracellular Matrix",
                                 '18' = "Astrocytes",
                                 '19' = "Unknown",
                                 '20' = "Deep Layer Neurons",
                                 '21' = "Erythroblasts")

DimPlot(CJ19wkCortexRep2)

#Save the names of the clusters
CJ19wkCortexRep2$celltype <- Idents(CJ19wkCortexRep2)

#Save this file so that it can be used for further analysis 
saveRDS(CJ19wkCortexRep2, file = '~/MarmosetWork2022/ScriptsMarmosetWork2022/CJ19wkRep2CortexCelltype.rds')

#Make custom color UMAP
CJ19wkUMAP <- DimPlot(CJ19wkCortexRep2, reduction = "umap", cols = c("Radial Glial Cells" = "magenta",
                                                             "Oligodendrocyte Progenitor Cells" = "chartreuse",
                                                             "Developing Oligodendrocytes" = "blue",
                                                             "Astrocytes" = "purple3",
                                                             "Unknown" = "gray",
                                                             "Deep Layer Neurons" = "orangered", 
                                                             "Ventral Intermediate Progenitors" = "springgreen4", 
                                                             "Endothelial Cells" = "black",
                                                             "Microglia" = "cyan", 
                                                             "Pericytes" = "olivedrab4",
                                                             "Dividing Progenitors" = "deepskyblue", 
                                                             "Extracellular Matrix" = "gold1",
                                                             "Erythroblasts"= "red4"))
CJ19wkUMAP
ggsave("~/MarmosetWork2022/CJCortex19wkRep2DataAnalysis/CJ19wkCortexRep2UMAPCustomColors.jpeg", width = 8, height = 4, units = c("in"), dpi = 300)