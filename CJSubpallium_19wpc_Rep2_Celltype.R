library(Seurat)
library(dplyr)
library(ggplot2)

#Load the previously saved Seurat file for CJ Subpallium
CJSubpallium <- readRDS(file = "~/MarmosetWork2022/ScriptsMarmosetWork2022/CJSubpallium19wkRep2Seur.rds")

#Identify the clusters based on the gene expression profiles using FindAllMarkers function of Seurat
#Renaming the clusters after marker analysis for each cluster

CJSubpallium <- RenameIdents(CJSubpallium,
                             '0' = "Radial Glial Cells", 
                             '1' = "Oligodendrocyte Progenitor Cells", 
                             '2' = "Mature Neurons",
                             '3' = "Immature Interneurons", 
                             '4' = "Radial Glial Cells", 
                             '5' = "Immature Interneurons", 
                             '6' = "Immature SST Interneurons", 
                             '7' = "Activated Microglia", 
                             '8' = "Endothelial Cells", 
                             '9' = "Mature Neurons",
                             '10' = "Pericytes", 
                             '11' = "Cycling Progenitors", 
                             '12' = "Activated Microglia", 
                             '13' = "Unknown", 
                             '14' = "Oligodendrocyte Progenitor Cells",
                             '15' = "Immature SST Interneurons",
                             '16' = "Microglia",
                             '17' = "Other",
                             '18' = "Erythroblasts",
                             '19' = "Mature Neurons",
                             '20' = "Ependymal Cells",
                             '21' = "Cajal-Retzius Neurons",
                             '22' = "Unknown")

#Save the names of the clusters
CJSubpallium$celltype <- Idents(CJSubpallium)

DimPlot(CJSubpallium)

#Save this file so that it can be used for further analysis 
saveRDS(CJSubpallium, file = '~/MarmosetWork2022/ScriptsMarmosetWork2022/CJ19wkRep2SubpalliumCelltype.rds')

#Make custom color UMAP
CJSubpallium19wkUMAP <- DimPlot(CJSubpallium, reduction = "umap",cols = c("Radial Glial Cells" = "magenta",
                                                                           "Oligodendrocyte Progenitor Cells" = "chartreuse",
                                                                           "Mature Neurons" = "gold1",
                                                                           "Immature Interneurons" = "orangered",
                                                                           "Unknown" = "gray",
                                                                           "Activated Microglia" = "purple3",
                                                                           "Endothelial Cells" = "black",
                                                                           "Microglia" = "cyan", 
                                                                           "Pericytes" = "olivedrab4",
                                                                           "Cycling Progenitors" = "deepskyblue", 
                                                                           "Immature SST Interneurons" = "blue",
                                                                           "Erythroblasts"= "red4",
                                                                           "Ependymal Cells" = "springgreen4", 
                                                                           "Cajal-Retzius Neurons" = "orchid", 
                                                                           "Other"= "gray48"))
CJSubpallium19wkUMAP
ggsave("~/MarmosetWork2022/CJSubpallium19wkRep2DataAnalysis/CJ19wkSubpalliumRep2UMAPCustomColors.jpeg", width = 8, height = 4, units = c("in"), dpi = 300)
