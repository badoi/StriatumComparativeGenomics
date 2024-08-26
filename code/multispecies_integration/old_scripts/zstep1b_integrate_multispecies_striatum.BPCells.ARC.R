library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(BPCells)
library(data.table)
library(future)
library(viridis)
library(harmony)
library(Nebulosa)
library(here)

# Enable parallelization
plan("sequential")
options(future.globals.maxSize = 100 * 1024^3, future.seed=TRUE)

DATADIR = 'data/tidy_data'
TEMPDIR = '/scratch/bnphan'
PLOTDIR = 'figures/multispecies_integration/plots'
dir.create(PLOTDIR, recursive = T)

############################
## 1) load in the datasets

# the names of obj_list in a vector written out by hand
list_names = c('Mouse_Stanley', # 1
               'Mouse_Saunders', # 2
               'Mouse_Chen', # 3
               'Mouse_Zeng', # 4
               'Rat_Phillips', # 5
               'Macaque_HeKleyman1', # 6
               'Macaque_HeKleyman2', # 7
               'Macaque_Chiou', # 8
               'Human_Tran', # 9
               'Human_Gayden', # 10
               'Human_Phan', # 11
               'Human_Siletti', # 12
               'Mouse_Zu', # 13
               'Human_Corces', # 14
               'Human_Li') # 15

# using the dataset orders to build the tree to integrate together
sample.tree = matrix(c(-1, -2, # 1: the mice, older data
                       -3, -4, # 2: the mice, newer data
                        1,  2, # 3: all mice RNA
                       -5,  3, # 4: all rodents RNA
                      -13,  4, # 5: all rodents RNA + ATAC
                       -6, -7, # 6: the monkeys, older RNA
                       -8,  6, # 7: all monkeys, RNA
                       -9,-10, # 8: the human, NAc
                      -11,  8, # 9: the human, NAc + dSTR
                      -12,  9, # 10: all human RNA
                        7, 10, # 11: all primates, RNA
                      -14,-15, # 12: all human ATAC
                       12, 11, # 13: all primates RNA + ATAC
                       13,  5  # 14: all mammals RNA + ATAC
), ncol = 2, byrow = T)

## copy the files to the temp directory for faster IO
files_list = list_names %>% paste0('.LogNorm.rds') %>% here(DATADIR, 'rdas/tmp', .)
files_list2 = list_names %>% paste0('.LogNorm.rds') %>% here(TEMPDIR, .)
thecall = paste('rsync -Paq', files_list, files_list2)
dir.create(TEMPDIR, recursive = T)
parallel::mclapply(thecall, system, mc.cores = 15)

## read in the sc transformed objects from the temp directory
obj_list <- parallel::mclapply(setNames(files_list2, list_names), readRDS, mc.cores = 15)

#######################################################
## 2) create BPCells counts directory for each dataset
BPCells_dir = list_names %>% paste0('_counts') %>% here(DATADIR, 'rdas/bpcells', .)
metadata_fn = list_names %>%  paste0('.metadata.rds') %>% here(DATADIR, 'rdas', .)
names(BPCells_dir) = names(metadata_fn) = list_names

lapply(list_names, function(x) {
  mat = convert_matrix_type(obj_list[[x]][["RNA"]]$counts, type ="uint32_t")
  write_matrix_dir(mat = mat, dir = BPCells_dir[x], overwrite = T)
})
lapply(list_names, function(x) {
  saveRDS(obj_list[[x]][[]], metadata_fn[x])
})

# now load in the counts and metadata lists
metadata.list = lapply(metadata_fn, readRDS) %>% bind_rows(.id = 'dataset')
data.list = lapply(BPCells_dir, open_matrix_dir)

cols_keep = c("Project", "Replicate", 'Assay', 'Species', 'Sex', 'Case', 
              'Region', 'Publication.Celltype', 'nCount_SCT', 'nFeature_SCT')
metadata <- lapply(metadata.list, function(x) x[, cols_keep, drop = F]) %>% 
  bind_rows(.id = 'dataset')

obj <- CreateSeuratObject(counts = data.list, meta.data = metadata)
obj

save_bp = here(DATADIR, 'rdas/merged_multispecies_striatum.bpcells.rds')
saveRDS(object = obj, file = save_bp)

# normalize the data, find variable features, and sketch the data per dataset
obj <- obj %>% NormalizeData() %>% FindVariableFeatures()
obj <- SketchData(obj, ncells = 5000, method = "LeverageScore", 
                  sketched.assay = "sketch")

## proceed with the integration
DefaultAssay(obj) <- "sketch"
obj <- obj %>% FindVariableFeatures(verbose = F) %>% 
  ScaleData(verbose = F) %>% RunPCA(verbose = F)

save_bp = here(DATADIR, 'rdas/merged_multispecies_striatum.bpcells.rds')
saveRDS(object = obj, file = save_bp)

# integrate the datasets at the sketch level
obj <- IntegrateLayers(obj, method = RPCAIntegration, orig = "pca", 
                       new.reduction = "integrated.rpca", verbose = T,
                       dims = 1:20, k.anchor = 25, sample.tree = sample.tree)

# cluster the sketched, integrated data
obj <- obj %>% 
  RunUMAP(reduction = "integrated.rpca", dims = 1:30, return.model = T) %>% 
  FindNeighbors(reduction = "integrated.rpca",dims = 1:30) %>% 
  FindClusters(resolution = 0.8, algorithm = 2)

## make the plots of this sketched, integrated dataset
pdf(here::here(PLOTDIR, 'Replicates_plots_bySpecies_clusters.sketched.umap.pdf'), 
    width = 11, height = 5)
DimPlot(obj, reduction = "umap", split.by = "Species", 
        cols = ArchR::paletteDiscrete(obj$seurat_clusters),
        label = T, raster = T,
        group.by = 'seurat_clusters') + theme(legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 2, override.aes = list(size = 3))) &
  FontSize(main = 7, x.title = 7, y.title = 7)
dev.off()

pdf(here::here(PLOTDIR, 'Replicates_plots_byProject_origLabel.sketched.umap.pdf'), 
    width = 30, height = 14)
DimPlot(obj, reduction = "umap", split.by = c("Project"), 
        cols = ArchR::paletteDiscrete(obj$Publication.Celltype),
        label = T, raster = F,
        group.by = 'Publication.Celltype') + theme(legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 3, override.aes = list(size = 3))) &
  FontSize(main = 4, x.title = 7, y.title = 7)
dev.off()

#############################################
##  project the sketched integration onto all cells
obj <- ProjectIntegration(object = obj, sketched.assay = "sketch", 
                          assay = "RNA", reduction = "integrated.rpca")

obj <- ProjectData(object = obj, sketched.assay = "sketch", assay = "RNA", 
                      sketched.reduction = "integrated.rpca",
                      full.reduction = "integrated.rpca.full", dims = 1:30, 
                      refdata = list(seurat_clusters.full = "seurat_clusters"))

obj <- RunUMAP(obj, reduction = "integrated.rpca.full", dims = 1:30, return.model = T, 
               reduction.name = "umap.full", reduction.key = "UMAP_full_")

save_bp = here(DATADIR, 'rdas/merged_multispecies_striatum.bpcells.rds')
saveRDS(object = obj, file = save_bp)

## make the plots of this full, integrated dataset
pdf(here::here(PLOTDIR, 'Replicates_plots_bySpecies_seuratClusters.full.umap.pdf'), 
    width = 8, height = 5)
DimPlot(obj, reduction = "umap.full", split.by = "Species", raster = T,
        label = T, group.by = 'seurat_clusters.full') +
  theme(legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 2, override.aes = list(size = 3))) &
  FontSize(main = 7, x.title = 7, y.title = 7)
dev.off()

## manually update clusters to be the main and sub clusters
obj@meta.data = obj[[]] %>% 
  mutate(
    integrated_clusters = case_when(
      seurat_clusters.full %in% c(1, 2, 3, 7, 15) ~ 'D2',  
      seurat_clusters.full %in% c(0, 4, 5, 6, 12) ~ 'D1', 
      seurat_clusters.full %in% c(9, 11) ~ 'eSPN', 
      seurat_clusters.full %in% c(13) ~ 'IC', 
      TRUE ~ 'Drop'), 
    
    integrated_subclusters = case_when(
      seurat_clusters.full %in% c(1, 2, 3) ~ 'D2-Matrix',  
      seurat_clusters.full %in% c(7, 15) ~ 'D2-Striosome',  
      seurat_clusters.full %in% c(0, 4, 6, 12) ~ 'D1-Matrix', 
      seurat_clusters.full %in% c(5) ~ 'D1-Striosome', 
      seurat_clusters.full %in% c(9) ~ 'D1/D2H', 
      seurat_clusters.full %in% c(11) ~ 'D1-NUDAP', 
      seurat_clusters.full %in% c(13) ~ 'IC', 
      TRUE ~ 'Drop')
  )

## relabel some things for organisation
proj_lvl = c('Siletti et al.', 'Tran et al.', 'Phan et al.', 'Li et al.',
             'Corces et al.', 
             'HeKleyman et al.', 'Chiou et al.', 'Zeng et al.', 
             'Stanley et al.',  'Saunders et al.', 'Chen et al.', 'Zu et al.',  
             'Phillips et al.', 'Savell et al.')

obj$Project = factor(obj$Project, proj_lvl)

levels(obj@meta.data$Project) = levels(obj$Project) %>%
  str_replace_all('HeK', 'He, K')

rep_lvl = obj[[]] %>% arrange(Project, Replicate) %>% 
  pull(Replicate) %>% unique()
obj$Replicate = factor(obj$Replicate, rep_lvl)

## make original cell type levels
cell_lvl = obj[[]] %>% arrange(Project, Publication.Celltype) %>% 
  pull(Publication.Celltype) %>% unique()
obj$Publication.Celltype = 
  factor(obj$Publication.Celltype, cell_lvl)

table(is.na(obj$integrated_clusters), obj$Project)
table(is.na(obj$integrated_clusters), obj$seurat_clusters.full)
table(is.na(obj$integrated_subclusters), obj$seurat_clusters.full)

## make the plots
cols = setNames(c(RColorBrewer::brewer.pal(7,'Paired'), 'black'), 
                c('D1-Matrix', 'D1-Striosome', 'D1-NUDAP',  'D1/D2H', 
                  'D2-Matrix', 'D2-Striosome', 'IC', 'Drop'))

obj@meta.data$integrated_subclusters = 
  factor(obj$integrated_subclusters, names(cols))

## make the plots
cols1 = obj[[]] %>% filter(integrated_clusters != 'Drop') %>% 
  pull(integrated_clusters) %>%  ArchR::paletteDiscrete()
cols1 = c(cols1, 'Drop' = 'black')
pdf(here::here(PLOTDIR, 'Replicates_plots_bySpecies_mainCluster.full.umap.pdf'), 
    width = 8, height = 5)
DimPlot(obj, reduction = "umap.full", split.by = "Species", 
        cols = cols1,
        label = T, group.by = 'integrated_clusters', raster = T) + 
  theme(legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 3))) &
  FontSize(main = 7, x.title = 7, y.title = 7)
dev.off()

pdf(here::here(PLOTDIR, 'Replicates_plots_bySpecies_subCluster.full.umap.pdf'), 
    width = 8, height = 5)
DimPlot(obj, reduction = "umap.full", split.by = "Species", raster = T,
        cols = cols, label = T, group.by = 'integrated_subclusters') +
  theme(legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 2, override.aes = list(size = 3))) &
  FontSize(main = 7, x.title = 7, y.title = 7)
dev.off()

## save the object
save_bp = here(DATADIR, 'rdas/integrated_multispecies_striatum.bpcells.rds')
saveRDS(object = obj, file = save_bp)

## convert the BPcells to a Matrix sparseMatrix
table(obj$integrated_subclusters, obj$integrated_clusters)

obj2 = subset(obj, subset = integrated_clusters != 'Drop') %>% JoinLayers()
DefaultAssay(obj2) = 'RNA' 


