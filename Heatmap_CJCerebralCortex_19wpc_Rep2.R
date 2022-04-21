library(Seurat)
library(dplyr)
library(ggplot2)

#Making heatmap for CJ Cortex 19wpc Rep2
CJCortex <- readRDS(file = '~/MarmosetWork2022/ScriptsMarmosetWork2022/CJ19wkRep2CortexCelltype.rds')

DimPlot(CJCortex)

#Find the levels of the cell types for the preferred order of cell types
y=levels(CJCortex)

#Preferred order of levels
#Radial Glial Cells, Cycling Progenitors, Ventral Intermediate Progenitors, Deep Layer Neurons, Oligodendrocyte Progenitor Cells, 
#Developing Oligodendrocytes, Astrocytes, Microglia, Pericytes, Endothelial Cells, 
#Extracellular Matrix, Erythroblasts, Unknown

mylevels1 = levels(CJCortex)[c(1,11,6,5,2,9,3,7,8,10,12,13,4 )]

Idents(CJCortex) <- factor(Idents(CJCortex), levels = mylevels1)

#Genes used in the heatmap
features <- c("HES1", "HES5", "NES", "EMX1", 
              "MKI67", "CENPE",
              "DLX1", "DLX5", "BCL11B", "ST18", 
              "PDGFRA", "OLIG1", "CLNDN11", 
              "AQP4", "S100B", "AIF1", "GAL",
              "RGS5", "COL4A1", "FLT1", "CLDN5", "APOD", "COL6A3", "ALAS2")

#Making heatmap
heatmapCJCortex <- DoHeatmap(CJCortex,
                          features = features,
                          assay = 'RNA',
                          group.by = "ident", 
                          slot = "data", 
                          lines.width = 3,
                          disp.min = 1.0,
                          disp.max = 2.0,
                          group.colors = c("magenta","deepskyblue","springgreen4","orangered","chartreuse","blue","purple3","cyan","olivedrab4","black","gold1","red4","gray" ),
                          group.bar.height = 0.05)
heatmapCJCortex + scale_fill_gradientn(colors = c("gray96",  "red"))

ggsave("~/MarmosetWork2022/CJCortex19wkRep2DataAnalysis/CJCortex19wkRep2Heatmap.tiff", width = 18, height = 8, units = c("in"), dpi = 300, bg= "white")
dev.off()

