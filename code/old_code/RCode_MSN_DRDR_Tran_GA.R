library(Seurat)
library(dplyr)
library(data.table)
library(ggplot2)

setwd("/Users/ghadareda/Desktop/UPitt-CNUP/Labs/Pfenning/Data/Human-Tran/Processed")
load("SCE_NAc-n8_tran-etal.rda")

pbmc <- as.Seurat(sce.nac.tran)
pbmc <- subset(pbmc, cellType %like% "MSN")

Idents(pbmc) <- "cellType"
DimPlot(pbmc, label = T)
DimPlot(pbmc,split.by = "sex" , label = T) + NoLegend()
VlnPlot(pbmc, features = "DRD3", pt.size = 0, split.by = "sex")
FeaturePlot(pbmc, "DRD3")


###############################

DefaultAssay(pbmc) <- "RNA"
VlnPlot(pbmc, "DRD1", slot = "data",pt.size = 0)
VlnPlot(pbmc, "DRD2", slot = "data", pt.size = 0)
VlnPlot(pbmc, "DRD3", slot = "data", pt.size = 0)

# Threshold 0.3 
Idents(pbmc) <- "NULL"

Idents(pbmc, WhichCells(object = pbmc, expression = DRD1 >= 1 & DRD2 >= 1 & DRD3 == 0, slot = 'counts')) <- "DRD1.DRD2"
Idents(pbmc, WhichCells(object = pbmc, expression = DRD1 >= 1 & DRD3 >= 1 & DRD2 == 0, slot = 'counts')) <- "DRD1.DRD3"
Idents(pbmc, WhichCells(object = pbmc, expression = DRD2 >= 1 & DRD3 >= 1 & DRD1 == 0, slot = 'counts')) <- "DRD2.DRD3"

Idents(pbmc, WhichCells(object = pbmc, expression = DRD1 >= 1 & DRD2 == 0 & DRD3 == 0, slot = 'counts')) <- "DRD1.only"
Idents(pbmc, WhichCells(object = pbmc, expression = DRD2 >= 1 & DRD1 == 0 & DRD3 == 0, slot = 'counts')) <- "DRD2.only"
Idents(pbmc, WhichCells(object = pbmc, expression = DRD3 >= 1 & DRD1 == 0 & DRD2 == 0, slot = 'counts')) <- "DRD3.only"

DRD1.DRD2 <- WhichCells(pbmc, idents = c("DRD1.DRD2"))
DRD1.DRD3 <- WhichCells(pbmc, idents = c("DRD1.DRD3"))
DRD2.DRD3 <- WhichCells(pbmc, idents = c("DRD2.DRD3"))

DRD1.only <- WhichCells(pbmc, idents = c("DRD1.only"))
DRD2.only <- WhichCells(pbmc, idents = c("DRD2.only"))
DRD3.only <- WhichCells(pbmc, idents = c("DRD3.only"))

Reduce(intersect, list(DRD1.DRD2,DRD1.DRD3,DRD2.DRD3))
Reduce(intersect, list(DRD1.only,DRD2.only,DRD3.only))
Reduce(intersect, list(DRD1.only,DRD2.only,DRD3.only,DRD1.DRD2,DRD1.DRD3,DRD2.DRD3))

## Plot
DimPlot(pbmc,cells.highlight= list(DRD1.DRD2,DRD2.DRD3, DRD1.DRD3), cols.highlight = c("darkblue", "darkred","darkgreen"), cols= "grey")
DimPlot(pbmc,cells.highlight= list(DRD1.only,DRD2.only,DRD3.only), cols.highlight = c("darkblue", "darkred","darkgreen"), cols= "grey")

DimPlot(pbmc,cells.highlight= list(DRD1.DRD2), cols.highlight = c("darkgreen"), cols= "grey")
DimPlot(pbmc,cells.highlight= list(DRD2.DRD3), cols.highlight = c("darkred"), cols= "grey")
DimPlot(pbmc,cells.highlight= list(DRD1.DRD3), cols.highlight = c("darkblue"), cols= "grey")

DimPlot(pbmc,cells.highlight= list(DRD1.only), cols.highlight = c("darkgreen"), cols= "grey")
DimPlot(pbmc,cells.highlight= list(DRD2.only), cols.highlight = c("darkred"), cols= "grey")
DimPlot(pbmc,cells.highlight= list(DRD3.only), cols.highlight = c("darkblue"), cols= "grey")

## Marker features
genes <- FindAllMarkers(pbmc, only.pos = T)
genes <- genes[genes$p_val_adj < 0.05,]
write.csv(genes,"/Users/ghadareda/Desktop/UPitt-CNUP/Labs/Pfenning/Projects/Striatum_Spatiomolecular_Map/Human-Tran/Data/Human_DRDR_markers.csv")

D1D2 <- genes[genes$cluster == "DRD1.DRD2",]
D1D3 <- genes[genes$cluster == "DRD1.DRD3",]
D2D3 <- genes[genes$cluster == "DRD2.DRD3",]
D1 <- genes[genes$cluster == "DRD1.only",]
D2 <- genes[genes$cluster == "DRD2.only",]
D3 <- genes[genes$cluster == "DRD3.only",]

# D1D2
D1D2_vs_D1D3 <- FindMarkers(pbmc, ident.1 = "DRD1.DRD2", ident.2 = "DRD1.DRD3", only.pos = T)
D1D2_vs_D2D3 <- FindMarkers(pbmc, ident.1 = "DRD1.DRD2", ident.2 = "DRD2.DRD3", only.pos = T)
D1D2_vs_D1D3 <- D1D2_vs_D1D3[genes$p_val_adj < 0.05,]
D1D2_vs_D2D3  <- D1D2_vs_D2D3[genes$p_val_adj < 0.05,]

# D1D3
D1D3_vs_D1D2 <- FindMarkers(pbmc, ident.1 = "DRD1.DRD3", ident.2 = "DRD1.DRD2", only.pos = T)
D1D3_vs_D2D3 <- FindMarkers(pbmc, ident.1 = "DRD1.DRD3", ident.2 = "DRD2.DRD3", only.pos = T)
D1D3_vs_D1D2 <- D1D3_vs_D1D2[genes$p_val_adj < 0.05,]
D1D3_vs_D2D3  <- D1D3_vs_D2D3[genes$p_val_adj < 0.05,]

# D2D3
D2D3_vs_D1D2 <- FindMarkers(pbmc, ident.1 = "DRD2.DRD3", ident.2 = "DRD1.DRD2", only.pos = T)
D2D3_vs_D1D3 <- FindMarkers(pbmc, ident.1 = "DRD2.DRD3", ident.2 = "DRD1.DRD3", only.pos = T)
D2D3_vs_D1D2 <- D2D3_vs_D1D2[genes$p_val_adj < 0.05,]
D2D3_vs_D1D3  <- D2D3_vs_D1D3[genes$p_val_adj < 0.05,]

# D1
D1_vs_D2 <- FindMarkers(pbmc, ident.1 = "DRD1.only", ident.2 = "DRD2.only", only.pos = T)
D1_vs_D3 <- FindMarkers(pbmc, ident.1 = "DRD1.only", ident.2 = "DRD3.only", only.pos = T)
D1_vs_D2 <- D1_vs_D2[genes$p_val_adj < 0.05,]
D1_vs_D3  <- D1_vs_D3[genes$p_val_adj < 0.05,]

# D2
D2_vs_D1 <- FindMarkers(pbmc, ident.1 = "DRD2.only", ident.2 = "DRD1.only", only.pos = T)
D2_vs_D3 <- FindMarkers(pbmc, ident.1 = "DRD2.only", ident.2 = "DRD3.only", only.pos = T)
D2_vs_D1 <- D2_vs_D1[genes$p_val_adj < 0.05,]
D2_vs_D3  <- D2_vs_D3[genes$p_val_adj < 0.05,]

# D3
D3_vs_D1 <- FindMarkers(pbmc, ident.1 = "DRD3.only", ident.2 = "DRD1.only", only.pos = T)
D3_vs_D2 <- FindMarkers(pbmc, ident.1 = "DRD3.only", ident.2 = "DRD2.only", only.pos = T)
D3_vs_D1 <- D3_vs_D1[genes$p_val_adj < 0.05,]
D3_vs_D2  <- D3_vs_D2[genes$p_val_adj < 0.05,]

# Overlapping genes 
D1D2_h <- intersect(intersect(D1D2$gene,rownames(D1D2_vs_D1D3)),rownames(D1D2_vs_D2D3))
D1D3_h <- intersect(intersect(D1D3$gene,rownames(D1D3_vs_D1D2)),rownames(D1D3_vs_D2D3))
D2D3_h <- intersect(intersect(D2D3$gene,rownames(D2D3_vs_D1D2)),rownames(D2D3_vs_D1D3))

D1_h <- intersect(intersect(D1$gene,rownames(D1_vs_D2)),rownames(D1_vs_D3))
D2_h <- intersect(intersect(D2$gene,rownames(D2_vs_D1)),rownames(D2_vs_D3))
D3_h <- intersect(intersect(D3$gene,rownames(D3_vs_D1)),rownames(D3_vs_D2))

# Conserved features in human and macaque
D1D2_hm <- intersect(D1D2_h,D1D2_m)
D1D3_hm <- intersect(D1D3_h,D1D3_m)
D2D3_hm <- intersect(D2D3_h,D2D3_m)

D1_hm <- intersect(D1_h,D1_m)
D2_hm <- intersect(D2_h,D2_m)
D3_hm <- intersect(D3_h,D3_m)

genes.toplot <- c("DRD1", "DRD2", "DRD3")
genes.toplot <- c("CHRM3" , "GRIK3" , "HS3ST4" ,"MSI2"  , "PTPRM" , "TTC12",
                  "CASZ1"   , "GRIK1"  ,  "GRIN2A" ,  "IVNS1ABP" ,  "OPCML", "OXR1"  ,   "SEMA6D" ,  "SORCS2" ,  "TACR1")
genes.toplot <- c("CHRM3" , "GRIK3" , "HS3ST4" ,"MSI2"  , "PTPRM" , "TTC12",
                  "FIG4", "PENK")
genes.toplot <- c("CNTN3", "EBF1", "RELN") 

genes.toplot <- c("DRD1", "DRD2", "DRD3",
                  "CHRM3" , "GRIK3" , "HS3ST4" ,"MSI2"  , "PTPRM" , "TTC12",
                  "CASZ1"   , "GRIK1"  ,  "GRIN2A" ,  "IVNS1ABP" ,  "OPCML", "OXR1"  ,   "SEMA6D" ,  "SORCS2" ,  "TACR1",
                  "FIG4", "PENK",
                  "CNTN3", "EBF1", "RELN") 

#pbmc@active.ident <- factor(x = pbmc@active.ident, levels = c("DRD1.only","DRD2.only","DRD3.only","DRD1.DRD2","DRD2.DRD3","DRD1.DRD3","NULL"))
#Idents(pbmc) <- factor(pbmc@active.ident, sort(levels(pbmc@active.ident)))
levels(pbmc) <- c("DRD1.DRD2","DRD2.DRD3","DRD1.DRD3","DRD1.only","DRD2.only","DRD3.only","NULL")

DotPlot(pbmc,features = genes.toplot, idents = c("DRD1.only","DRD2.only","DRD1.DRD2"), scale = T, cols = c("green", "blue")) + coord_flip()
DotPlot(pbmc,features = genes.toplot, idents = c("DRD2.only","DRD3.only","DRD2.DRD3"), scale = T, cols = c("green", "blue")) + coord_flip()
DotPlot(pbmc,features = genes.toplot, idents = c("DRD1.only","DRD3.only","DRD1.DRD3"),scale = F, cols = c("green", "blue")) + coord_flip()


VlnPlot(pbmc, "DRD1",idents = c("DRD1.DRD2","DRD2.DRD3","DRD1.DRD3","DRD1.only","DRD2.only","DRD3.only"),pt.size = 0)

##############################
