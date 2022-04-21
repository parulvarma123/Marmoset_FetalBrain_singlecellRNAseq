library(Seurat)
library(dplyr)
library(ggplot2)

#Making heatmap for CJ Subpallium 19wpc Rep2
CJSubpallium <- readRDS(file = '~/MarmosetWork2022/ScriptsMarmosetWork2022/CJ19wkRep2SubpalliumCelltype.rds')

DimPlot(CJSubpallium)

#Find the levels of the cell types for the preferred order of cell types
y=levels(CJSubpallium)

#Preferred order of levels
#levels = Radial Glial cells, Cycling Progenitors, Cajal-Retzius Neurons, Immature Interneurons, Immature SST Interneurons
#Mature neurons, Oligodendrocyte progenitor cells, microglia, activated microglia, pericytes, endothelial cells, 
#ependymal cells, erythroblasts, other, unknown

mylevels1 = levels(CJSubpallium)[c(1,9,15,4,5,3,2, 11, 6, 8, 7, 14, 13,12, 10 )]

Idents(CJSubpallium) <- factor(Idents(CJSubpallium), levels = mylevels1)

#Genes used in the heatmap
features <- c("HES1", "HES5", "SLC1A2", "AQP4", "GFAP", 
              "MKI67", "BIRC5", "LHX9",
              "MEF2C", 
              "DLX1", "DLX2", "DLX5","GAD1", "GAD2", "SP9", "NPAS1", "ARX", "MAF",
              "MAFB", "SST",
              "NEUROD6", "NRGN", 
              "PDGFRA", "OLIG1", 
              "TREM2", "CCL24", 
              "RGS5", 
              "CLDN5", "MSFD2A", "CXCR4",
              "HOATZ", "FAM183A",
              "ALAS2", "CD52")

#Making heatmap
heatmapCJSub <- DoHeatmap(CJSubpallium,
                          features = features,
                          assay = 'RNA',
                          group.by = "ident", 
                          slot = "data", 
                          lines.width = 3,
                          disp.min = 1.0,
                          disp.max = 2.0,
                          group.colors = c("magenta","deepskyblue","orchid", "orangered", "blue", "gold1","chartreuse","cyan","purple3","olivedrab4","black","springgreen4","red4","gray48","gray"),
                          group.bar.height = 0.05)
heatmapCJSub + scale_fill_gradientn(colors = c("gray96",  "red"))

ggsave("~/MarmosetWork2022/CJSubpallium19wkRep2DataAnalysis/CJSubpallium19wkRep2Heatmap.tiff", width = 20, height = 8, units = c("in"), dpi = 300, bg= "white")
dev.off()

