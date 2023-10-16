# load libraries
library(Seurat)
#library(SeuratDisk)
library(dplyr)
library(ggplot2)

RDA_PATH = "/projects/pfenninggroup/singleCell/Tran2021_human_snRNA-seq/SCE_NAc-n8_tran-etal.rda"
pbmc = load(RDA_PATH) 

if ("umap" %in% names(pbmc)) {
  # UMAP dimensionality reduction exists
  umap_coords <- pbmc$umap
  print(umap_coords)  # Print or perform any desired operation with the UMAP coordinates
} else {
  # UMAP dimensionality reduction does not exist
  print("UMAP dimensionality reduction not found.")
}

# get SingleCellExperiment Object
sce <- as(sce.nac.tran, "SingleCellExperiment")
sce
# Convert SingleCellExperiment to Seurat object
seurat_obj <- CreateSeuratObject(counts = assays(sce)$counts)

# Get cell type information from colData
cell_types <- colData(sce)$cellType

# Filter cell types that begin with "MSN"
filtered_cell_types <- cell_types[grep("^MSN", cell_types)]

# Filter the SingleCellExperiment object based on cell types
filtered_sce <- sce[, cell_types %in% filtered_cell_types]

# Reorder cell names in Seurat object to match the order in filtered_sce
seurat_obj <- seurat_obj[, colnames(seurat_obj) %in% colnames(filtered_sce)]

# Add cell type information to Seurat object
seurat_obj$cellType <- filtered_sce$cellType
# make plots
violin_plot <- VlnPlot(seurat_obj, features="DRD1", group.by = "cellType")
ggsave("Human_NAc_DRD1.png", plot= violin_plot, dpi= 300)
violin_plot <- VlnPlot(seurat_obj, features="DRD2", group.by = "cellType")
ggsave("Human_NAc_DRD2.png", plot= violin_plot, dpi= 300)
violin_plot <- VlnPlot(seurat_obj, features="DRD3", group.by = "cellType")
ggsave("Human_NAc_DRD3.png", plot= violin_plot, dpi= 300)

# Threshold 0.3 
Idents(seurat_obj, WhichCells(object = seurat_obj, expression = DRD1 > 0.3, slot = 'data')) <- "DRD1.pos"
Idents(seurat_obj, WhichCells(object = seurat_obj, expression = DRD2 > 0.3, slot = 'data')) <- "DRD2.pos"
Idents(seurat_obj, WhichCells(object = seurat_obj, expression = DRD3 > 0.3, slot = 'data')) <- "DRD3.pos"

Idents(seurat_obj, WhichCells(object = seurat_obj, expression = DRD1 < 0.3, slot = 'data')) <- "DRD1.neg"
Idents(seurat_obj, WhichCells(object = seurat_obj, expression = DRD2 < 0.3, slot = 'data')) <- "DRD2.neg"
Idents(seurat_obj, WhichCells(object = seurat_obj, expression = DRD3 < 0.3, slot = 'data')) <- "DRD3.neg"

Idents(seurat_obj, WhichCells(object = seurat_obj, expression = DRD1 > 0.3 & DRD2 > 0.3, slot = 'data')) <- "DRD1.DRD2.pos"
Idents(seurat_obj, WhichCells(object = seurat_obj, expression = DRD1 > 0.3 & DRD3 > 0.3, slot = 'data')) <- "DRD1.DRD3.pos"
Idents(seurat_obj, WhichCells(object = seurat_obj, expression = DRD2 > 0.3 & DRD3 > 0.3, slot = 'data')) <- "DRD2.DRD3.pos"

DRD1.DRD2 <- WhichCells(seurat_obj, idents = c("DRD1.DRD2.pos"))
DRD1.DRD3 <- WhichCells(seurat_obj, idents = c("DRD1.DRD3.pos"))
DRD2.DRD3 <- WhichCells(seurat_obj, idents = c("DRD2.DRD3.pos"))


reduction_method <- "umap"  # or "tsne" or "pca"

# Compute the dimensionality reduction
seurat_obj <- RunUMAP(seurat_obj, dimred="PCAopt")  # or RunTSNE(seurat_obj) or RunPCA(seurat_obj)

# Create a named list of cell types and their corresponding cell names
highlight_cells <- list(
  "DRD1.DRD2.pos" = WhichCells(seurat_obj, idents = "DRD1.DRD2.pos"),
  "DRD1.DRD3.pos" = WhichCells(seurat_obj, idents = "DRD1.DRD3.pos"),
  "DRD2.DRD3.pos" = WhichCells(seurat_obj, idents = "DRD2.DRD3.pos")
)

# Create a DimPlot with reduction method and highlight specific cell types
DimPlot(seurat_obj, reduction_method, cells.highlight = highlight_cells,
        cols.highlight = c("darkblue", "darkred", "darkgreen"), cols = "grey")


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


