# load libraries
library(Seurat)
#library(SeuratDisk)
library(dplyr)
library(ggplot2)

RDS_PATH = "/projects/pfenninggroup/singleCell/Savell2020_rat_snRNA-seq/data/tidy_data/rdas/Savell2020_rawCounts_snRNA-seq_N4.sce.rds" # fill Savell project path in here
pbmc = readRDS(RDS_PATH) # 2 monkeys / 2000 genes / 7000 cells
# how many rats, genes, cells?
# update object
sce <- as(pbmc, "SingleCellExperiment")
sce
seurat_obj <- CreateSeuratObject(counts = assays(sce)$counts)
cell_types <- colData(sce)$CellType

# Filter cell types that begin with "MSN"
filtered_cell_types <- cell_types[grep("MSN", cell_types)]

# Filter the SingleCellExperiment object based on cell types
filtered_sce <- sce[, cell_types %in% filtered_cell_types]

# Reorder cell names in Seurat object to match the order in filtered_sce
seurat_obj <- seurat_obj[, colnames(seurat_obj) %in% colnames(filtered_sce)]

# Add cell type information to Seurat object
seurat_obj$cellType <- filtered_sce$CellType

violin_plot <- VlnPlot(seurat_obj, features="Ddr1", group.by = "cellType")
ggsave("Rat_NAc_DRD1.png", plot= violin_plot, dpi= 300)
violin_plot <- VlnPlot(seurat_obj, features="Drd2", group.by = "cellType")
ggsave("Rat_NAc_DRD2.png", plot= violin_plot, dpi= 300)
violin_plot <- VlnPlot(seurat_obj, features="Drd3", group.by = "cellType")
ggsave("Rat_NAc_DRD3.png", plot= violin_plot, dpi= 300)

violin_plot <- VlnPlot(seurat_obj, features="Drd3", group.by = "Sex")
ggsave("Rat_NAc_DRD3_sex.png", plot= violin_plot, dpi= 300)

DefaultAssay(pbmc) <- "RNA"
Idents(pbmc) <- "MSN_type"

DefaultAssay(pbmc) <- "integrated"
DimPlot(pbmc)

DefaultAssay(pbmc) <- "RNA"
VlnPlot(pbmc, "DRD1", slot = "data")
VlnPlot(pbmc, "DRD2", slot = "data")
VlnPlot(pbmc, "DRD3", slot = "data")

# Threshold 0.3 
# where did this threshold come from? 
DefaultAssay(pbmc) <- "RNA"
Idents(seurat_obj, WhichCells(object = seurat_obj, expression = Drd1 > 0.3, slot = 'data')) <- "DRD1.pos"
Idents(seurat_obj, WhichCells(object = seurat_obj, expression = Drd2 > 0.3, slot = 'data')) <- "DRD2.pos"
Idents(seurat_obj, WhichCells(object = seurat_obj, expression = Drd3 > 0.3, slot = 'data')) <- "DRD3.pos"

Idents(seurat_obj, WhichCells(object = seurat_obj, expression = Drd1 < 0.3, slot = 'data')) <- "DRD1.neg"
Idents(seurat_obj, WhichCells(object = seurat_obj, expression = Drd2 < 0.3, slot = 'data')) <- "DRD2.neg"
Idents(seurat_obj, WhichCells(object = seurat_obj, expression = Drd3 < 0.3, slot = 'data')) <- "DRD3.neg"

Idents(seurat_obj, WhichCells(object = seurat_obj, expression = Drd1 > 0.3 & Drd2 > 0.3, slot = 'data')) <- "DRD1.DRD2.pos"
Idents(seurat_obj, WhichCells(object = seurat_obj, expression = Drd1 > 0.3 & Drd3 > 0.3, slot = 'data')) <- "DRD1.DRD3.pos"
Idents(seurat_obj, WhichCells(object = seurat_obj, expression = Drd2 > 0.3 & Drd3 > 0.3, slot = 'data')) <- "DRD2.DRD3.pos"

DRD1.DRD2 <- WhichCells(seurat_obj, idents = c("DRD1.DRD2.pos"))
DRD1.DRD3 <- WhichCells(seurat_obj, idents = c("DRD1.DRD3.pos"))
DRD2.DRD3 <- WhichCells(seurat_obj, idents = c("DRD2.DRD3.pos"))

DimPlot(seurat_obj, "umap", cells.highlight= list(DRD1.DRD2,DRD1.DRD3,DRD2.DRD3), cols.highlight = c("darkblue", "darkred","darkgreen"), cols= "grey")


genes <- FindAllMarkers(pbmc, only.pos = T)
write.csv(genes,"/Users/ghadareda/Desktop/UPitt-CNUP/Labs/Pfenning/Data/DRDR_markers.csv")

genes2 <- FindMarkers(pbmc, ident.1 = "DRD2.pos", ident.2 = "DRD1.neg" , only.pos = T)
write.csv(genes2,"/Users/ghadareda/Desktop/UPitt-CNUP/Labs/Pfenning/Data/DRD2pos_DRD1neg_markers.csv")

VlnPlot(pbmc, "DRD1")
VlnPlot(pbmc, "DRD2")
VlnPlot(pbmc, "DRD3")

metadata <- pbmc@meta.data
metadata$ID_barcode <- rownames(metadata)
metadat <- metadata[,"ID_barcode",drop = FALSE]

DRDR_type <- as.data.frame(pbmc@active.ident)
DRDR_type$ID_barcode <- rownames(DRDR_type)

# D2 MSNs - GWAS


