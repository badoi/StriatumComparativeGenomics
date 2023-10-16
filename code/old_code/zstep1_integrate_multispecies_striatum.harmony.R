library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(harmony)
library(data.table)
library(future)
library(viridis)
library(Nebulosa)
library(here)

# Enable parallelization
plan("multicore", workers = 24)
options(future.globals.maxSize = 100 * 1024^3, future.seed=TRUE)

DATADIR = 'data/tidy_data/'
PLOTDIR = 'figures/multispecies_integration/plots'
dir.create(PLOTDIR, recursive = T)

##############################
## 1) load in the datasets

# mouse dataset w/ human genes
mouse_fn = here(DATADIR,'Mouse_Stanley/Stanley2020_mouse_striatum_msns.hgGenes.h5Seurat')
obj_mm = LoadH5Seurat(mouse_fn)
table(obj_mm$continuous.subtype, obj_mm$discrete.type)

obj_mm@meta.data = obj_mm[[]] %>% 
  mutate(Region = case_when(
    grepl('^Pat|^Mat|CPu|pat', continuous.subtype) ~ 'CPu',
    grepl('ICj|^OT|vClust|unknown', continuous.subtype) ~ 'OT',
    grepl('Shell|NAcc', continuous.subtype) ~ 'NAcc'), 
    Species = "Mouse", Project = 'Stanley_et_al',
    celltype.orig = discrete.type, 
    Replicate = expt.name)

table(obj_mm$continuous.subtype, obj_mm$Region)

# Tran et al. NAcc human dataset
tran_fn = here(DATADIR,'Human_Tran/SCE_NAc-n8_tran-etal.msns.h5Seurat')
obj_hgNAc = LoadH5Seurat(tran_fn)
obj_hgNAc = RenameAssays(object = obj_hgNAc, originalexp = 'RNA')
obj_hgNAc$Region = "NAcc"
obj_hgNAc$Species = "Human"
obj_hgNAc$Project = "Tran_et_al"
obj_hgNAc$Replicate = obj_hgNAc$donor
obj_hgNAc$celltype.orig = as.character(obj_hgNAc$cellType)

# Phan et al. Caud/Put human CTL only dataset
phan_fn = here(DATADIR,'Human_Phan/OUD_Striatum_refined_SeuratObj_N10.msns.h5Seurat')
obj_hgDStr = LoadH5Seurat(phan_fn)
DefaultAssay(obj_hgDStr) = 'RNA'
obj_hgDStr$Species = "Human"
obj_hgDStr$Project = "Phan_et_al"
obj_hgDStr$Replicate = obj_hgDStr$ID
obj_hgDStr$celltype.orig = obj_hgDStr$celltype3

# He, Kleyman et al. Caud/Put/Nacc macaque dataset
macaque_fn = here(DATADIR,'Macaque_He/GSE167920_Results_MSNs_processed_final.h5Seurat')
obj_rm = LoadH5Seurat(macaque_fn)
DefaultAssay(obj_rm) = 'RNA'
obj_rm@meta.data = obj_rm[[]] %>% 
  mutate(Region = case_when(
    region_name == 'caudate' ~ "Caudate", 
    region_name == 'putamen' ~ 'Putamen', 
    region_name == 'nacc' ~ 'NAcc'), 
    Species = "Macaque", Project = "HeKleyman_et_al",
    Replicate = monkey,
    celltype.orig = as.character(MSN_type))

####################################################################
## 2) integrate the datasets by shared common genes
obj_list = list('Human_Tran' = obj_hgNAc, 'Human_Phan' = obj_hgDStr, 
                'Macaque_HeKleyman' = obj_rm, 'Mouse_Stanley' = obj_mm)

# subset to just the genes that overlap all datasets
genes = obj_list %>% lapply(rownames) %>% Reduce(f = 'intersect')
length(genes) #11740 genes

## concatenate together all the datasets and perform SCTransform
obj_list = obj_list %>% lapply(function(x) x[genes,]) # subset to just shared genes
obj_merge <- merge(obj_list[[1]], y = obj_list[-1], 
                  add.cell.ids = names(obj_list), project = "Multispecies Striatum")
obj_merge = obj_merge %>% SCTransform(method = "glmGamPoi") %>% 
  RunPCA(verbose = F, assay = 'SCT') 
obj_merge = obj_merge %>%
  RunHarmony(c("Project", 'Species', 'Replicate'),max.iter.harmony = 20, 
             lambda = c(.2, .2, .2))  %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.1, algorithm = 1)

# ## manually update clusters
# obj_merge@meta.data = obj_merge[[]] %>% 
#   mutate(integrated_clusters = case_when(seurat_clusters == 0 ~ 'D2', 
#                                          seurat_clusters == 1 ~ 'D1', 
#                                          seurat_clusters == 2 ~ 'D1/2H', 
#                                          grepl('ICj', celltype.orig) ~ 'ICj', 
#                                          seurat_clusters == 4 ~ 'ICj', 
#                                          seurat_clusters == 5 ~ 'MSN_human', 
#                                          seurat_clusters == 6 ~ 'MSN_primA', 
#                                          seurat_clusters == 7 ~ 'MSN_primB', 
#   ))

## save the integrated object
save_fn = here(DATADIR,'rdas/integrated_mouse_macaque_human_striatum.harmony.h5Seurat')
SaveH5Seurat(obj_merge, save_fn, overwrite = T)


###################################
## 3) make QC visualization plots

# grab the marker genes that should be plot later on
markers =  c('DRD1','DRD2', 'DRD3', 'FOXP2', 'PDYN', 'PENK', 'OPRM1', 'OPRK1',
             'CPNE4', 'CASZ1', 'OTOF', 'CACNG5')
markers = markers[markers %in% genes]

cols = ArchR::paletteDiscrete(obj_merge$celltype.orig)
cols2 = ArchR::paletteDiscrete(obj_merge$seurat_clusters, set = 'stallion2')

pdf(here(PLOTDIR, 'original_clusters_byProject.pdf'), width = 8.5, height = 5)
DimPlot(obj_merge, reduction = "umap", split.by = "Project", cols = cols, 
        group.by = 'celltype.orig', label = T, label.size = 3) +
  guides(colour = guide_legend(nrow = 3, override.aes = list(size = 2))) +
  theme(legend.position = 'bottom', legend.text=element_text(size=7)) + 
  ggtitle("Multi-species integrated datasets: original cluster labels")
dev.off()


pdf(here(PLOTDIR, 'integrated_clusters_byProject.pdf'), width = 8.5, height = 5)
DimPlot(obj_merge, reduction = "umap", split.by = "Project", cols = cols2,
        group.by = 'seurat_clusters', label = T, label.size = 5) +
  theme(legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 4))) + 
  ggtitle("Multi-species integrated datasets: integrated cluster labels")
dev.off()


pdf(here(PLOTDIR, 'DRD_expression_byProject.pdf'), width = 8.5, height = 5)
FeaturePlot(obj_merge, c("DRD1", "DRD2", 'DRD3')) & 
  viridis::scale_color_viridis(option="magma")
dev.off()


pdf(here(PLOTDIR, 'QC_plots_byProject.pdf'), width = 8.5, height = 5)
FeaturePlot(obj_merge, c("nCount_SCT", "nFeature_SCT"), split.by = "Project") & 
  viridis::scale_color_viridis(option="magma")

DimPlot(obj_merge, reduction = "umap", split.by = "Project", 
        cols = ArchR::paletteDiscrete(obj_merge$Replicate),
        group.by = 'Replicate', label = T, label.size = 5)
dev.off()


