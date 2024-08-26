library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(data.table)
library(future)
library(viridis)
library(harmony)
library(Nebulosa)
library(here)

# Enable parallelization
plan("multicore", workers = 6)
options(future.globals.maxSize = 100 * 1024^3, future.seed=TRUE)

DATADIR = 'data/tidy_data/'
PLOTDIR = 'figures/multispecies_integration/plots'
dir.create(PLOTDIR, recursive = T)

threshold = 0.5

cols = ArchR::paletteDiscrete(c('D1', 'D2', 'eSPN', 'IC'), set = 'paired')

cols2 = setNames(RColorBrewer::brewer.pal(7,'Paired'), 
                c('D1-Matrix', 'D1-Striosome', 'D1-NUDAP',  'D1/D2H', 
                  'D2-Matrix', 'D2-Striosome', 'IC'))

##################################
## 1) lead the integrated RNA object
save_fn = here(DATADIR,'rdas/integrated_mouse_rat_macaque_human_striatum.RNA.rds')
combined.sct = LoadSeuratRds(save_fn)
DefaultAssay(combined.sct) = 'integrated'

##############################
## 2) load in the  Corces et al.
# Donor information from https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-020-00721-x/MediaObjects/41588_2020_721_MOESM3_ESM.xlsx
corces_fn = here(DATADIR,'Human_Corces/Corces2020_human_striatum_msns.hgGenes.h5Seurat')
obj_hg = LoadH5Seurat(corces_fn)
obj_hg$Species = "Human"
obj_hg$Project = "Corces_et_al"
obj_hg$Region = "Caudate"
obj_hg$Replicate = obj_hg$Subject
obj_hg$Publication.Celltype = obj_hg$ClusterName
obj_hg$Case = 'Unaffected'
obj_hg$Sex = ifelse(obj_hg$Replicate %in% c('09_1589', '14_1018'), 'Female', 'Male')

table(obj_hg$Sex, obj_hg$Replicate)


## grab the human cells from the integrated object
ref_human = combined.sct %>% subset(Project == 'Phan et al.')
DefaultAssay(ref_human) = 'SCT'

ref_human = ref_human %>% FindVariableFeatures(assay = 'SCT') %>% 
  RunPCA(verbose = F, assay = 'SCT') %>% 
  RunHarmony(c("orig.ident", 'Sex', 'Region', 'DSM.IV.OUD'), 
             max_iter = 20, lambda = rep(1,4)) %>% 
  RunUMAP(reduction = "harmony", dims = 1:30, return.model = T) 

## map the query objects to the reference
anchors <- FindTransferAnchors(
  reference = ref_human, query = obj_hg, normalization.method = "SCT", l2.norm = T,
  k.anchor = 40, dims = 1:30, reference.reduction = "pca")

obj_hg <- MapQuery(
  anchorset = anchors, reference = ref_human,  query = obj_hg,
  reference.reduction = "pca", reduction.model = "umap",
  refdata = list(integrated_clusters = "integrated_clusters",
                 integrated_subclusters = "integrated_subclusters" ))

## make the plots
pdf(here(PLOTDIR, 'Corces2020_human_striatum_msns.projected.pdf'), 
    width = 10, height = 5)

DimPlot(obj_hg, reduction = "ref.umap", label = T, label.size = 3, 
        split.by = "Replicate", repel = T, cols = cols,
        group.by = "predicted.integrated_clusters") + 
  NoLegend() + ggtitle("Query transferred integrated_clusters")

VlnPlot(obj_hg, features = c("DRD1", "DRD2", "DRD3", "TSHZ1", "PDYN", "PENK"), 
        group.by = "predicted.integrated_clusters", cols = cols, pt.size = 0)

DimPlot(obj_hg, reduction = "ref.umap", label = T, label.size = 3, 
        split.by = "Replicate", repel = T, cols = cols2,
        group.by = "predicted.integrated_subclusters") + 
  NoLegend() + ggtitle("Query transferred integrated_subclusters")

VlnPlot(obj_hg, features = c("DRD1", "DRD2", "DRD3", "TSHZ1", "PDYN", "PENK"), 
        group.by = "predicted.integrated_subclusters", cols = cols2, pt.size = 0)

dev.off()

corces_fn2 = here(DATADIR,'Human_Corces/Corces2020_human_striatum_msns.hgGenes.annotated.rds')
SaveSeuratRds(obj_hg, corces_fn2)










# 3, Zu et al.CATlas BICCN mouse
obj_mm = here(DATADIR,'Mouse_Zu/Zu2023_mouse_striatum_msns.hgGenes.rds') %>% 
  readRDS()

obj_mm@meta.data = obj_mm[[]] %>% 
  mutate(Region = case_when(
    grepl('CP', SubRegion) ~ 'CPu',
    grepl('ACB', SubRegion) ~ 'NAcc/OT'), 
    Species = "Mouse", 
    Project = 'Zu_et_al',
    Publication.Celltype = NA,
    Replicate = Sample, 
    Sex = 'Male',
    Case = 'Unaffected')

table(obj_mm$Publication.Celltype)

## grab the mouse cells from the integrated object
ref_mouse = combined.sct %>% 
  subset(Project %in%  c('Stanley et al.', 'Saunders et al.')) %>% 
  RunPCA(verbose = F)

## map the query objects to the reference
anchors <- FindTransferAnchors(
  reference = ref_mouse, query = obj_mm, normalization.method = "SCT", 
  l2.norm = T, k.anchor = 40, dims = 1:30, reference.reduction = "pca")

obj_mm <- MapQuery(
  anchorset = anchors, reference = ref_mouse,  query = obj_mm,
  reference.reduction = "pca", reduction.model = "umap",
  refdata = list(integrated_clusters = "integrated_clusters",
                 integrated_subclusters = "integrated_subclusters" ))

## make the plots
pdf(here(PLOTDIR, 'Zu2023_mouse_striatum_msns.projected.pdf'), 
    width = 10, height = 5)

DimPlot(obj_mm, reduction = "ref.umap", label = T, label.size = 3, 
        split.by = "Region", repel = T, cols = cols,
        group.by = "predicted.integrated_clusters") + 
  NoLegend() + ggtitle("Query transferred integrated_clusters")

VlnPlot(obj_mm, features = c("DRD1", "DRD2", "DRD3", "TSHZ1", "PDYN", "PENK"), 
        group.by = "predicted.integrated_clusters", cols = cols, pt.size = 0)

DimPlot(obj_mm, reduction = "ref.umap", label = T, label.size = 3, 
        split.by = "Region", repel = T, cols = cols2,
        group.by = "predicted.integrated_subclusters") + 
  NoLegend() + ggtitle("Query transferred integrated_subclusters")

VlnPlot(obj_mm, features = c("DRD1", "DRD2", "DRD3", "TSHZ1", "PDYN", "PENK"), 
        group.by = "predicted.integrated_subclusters", cols = cols2, pt.size = 0)
dev.off()


zu_fn2 = here(DATADIR,'Mouse_Zu/Zu2023_mouse_striatum_msns.hgGenes.annotated.rds')
SaveSeuratRds(obj_mm, zu_fn2)

