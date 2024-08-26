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
plan("sequential")
options(future.globals.maxSize = 100 * 1024^3, future.seed=TRUE)

DATADIR = 'data/tidy_data/'
TEMPDIR = '/scratch/bnphan'
PLOTDIR = 'figures/multispecies_integration/plots'
dir.create(PLOTDIR, recursive = T)

##############################
## 1) load in the datasets

# the names of obj_list in a vector written out by hand
list_names = c('Mouse_Stanley', 'Mouse_Saunders', 'Mouse_Chen', 'Mouse_Zeng', 
               'Rat_Phillips', 'Macaque_HeKleyman1', 'Macaque_HeKleyman2', 
               'Macaque_Chiou', 'Human_Phan', 'Human_Tran', 'Human_Siletti', 
               'Mouse_Zu', 'Human_Corces', 'Human_Li')

## copy the files to the temp directory for faster IO
files_list = list_names %>% paste0('.sctransformed.rds') %>% here(DATADIR, 'rdas/tmp', .)
files_list2 = list_names %>% paste0('.sctransformed.rds') %>% here(TEMPDIR, .)
thecall = paste('rsync -Paq', files_list, files_list2)
dir.create(TEMPDIR, recursive = T)
parallel::mclapply(thecall, system, mc.cores = 14)

## read in the sc transformed objects from the temp directory
obj_list <- parallel::mclapply(setNames(files_list2, list_names), readRDS, mc.cores = 14)

## proceed with the integration
features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = 3000)
obj_list <- PrepSCTIntegration(object.list = obj_list, anchor.features = features)
obj_list <- lapply(X = obj_list, FUN = RunPCA, features = features, verbose = F)
pryr::mem_used()

## integrate normalized objects
anchors <- FindIntegrationAnchors(obj_list, normalization.method = "SCT",
                                  reference = c(4, 11),
                                  reduction = "rpca", k.anchor = 40,
                                  anchor.features = features, dims = 1:30)

## do the integration of multiple datasets together
combined.sct = IntegrateData(anchorset = anchors, 
                             normalization.method = "SCT", 
                             dims = 1:30)

# use Harmony to finish out the batch correction for different samples
combined.sct = combined.sct %>% 
  RunPCA(verbose = F, assay = 'integrated') %>% 
  RunHarmony(c("Project", "Replicate", 'Assay', 
               'Species', 'Sex', 'Case', 'Region'), 
             max_iter = 20, lambda = rep(1,7)) %>% 
  RunUMAP(reduction = "harmony", dims = 1:30, return.model = T) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:30) %>% 
  # test this is the best resolution b/t lump & split
  FindClusters(resolution = 1, algorithm = 2)

###############################################################
## 3) plot the seurat clusters to see which clusters are which

## relabel some things for organisation
proj_lvl = c('Siletti et al.', 'Tran et al.', 'Phan et al.', 'Li et al.',
             'Corces et al.', 'HeKleyman et al.', 'Zeng et al.',
             'Stanley et al.',  'Saunders et al.', 'Chen et al.', 'Zu et al.',  
             'Phillips et al.', 'Savell et al.')
combined.sct$Project = factor(combined.sct$Project, proj_lvl)

levels(combined.sct@meta.data$Project) = levels(combined.sct$Project) %>%
  str_replace_all('HeK', 'He, K')

rep_lvl = combined.sct[[]] %>% arrange(Project, Replicate) %>% 
  pull(Replicate) %>% unique()
combined.sct$Replicate = factor(combined.sct$Replicate, rep_lvl)

## make original cell type levels
cell_lvl = combined.sct[[]] %>% arrange(Project, Publication.Celltype) %>% 
  pull(Publication.Celltype) %>% unique()
combined.sct$Publication.Celltype = 
  factor(combined.sct$Publication.Celltype, cell_lvl)

## make the plots
pdf(here::here(PLOTDIR, 'Replicates_plots_bySpecies_clusters.umap.pdf'), 
    width = 11, height = 5)
DimPlot(combined.sct, reduction = "umap", split.by = "Species", 
        cols = ArchR::paletteDiscrete(combined.sct$seurat_clusters),
        label = T, raster = F,
        group.by = 'seurat_clusters') + theme(legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 2, override.aes = list(size = 3))) &
  FontSize(main = 7, x.title = 7, y.title = 7)
dev.off()

pdf(here::here(PLOTDIR, 'Replicates_plots_byProject_origLabel.umap.pdf'), 
    width = 20, height = 7)
DimPlot(combined.sct, reduction = "umap", split.by = c("Project"), 
        cols = ArchR::paletteDiscrete(combined.sct$Publication.Celltype),
        label = T, raster = F,
        group.by = 'Publication.Celltype') + theme(legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 3, override.aes = list(size = 3))) &
  FontSize(main = 7, x.title = 7, y.title = 7)
dev.off()

## manually update clusters to be the main and sub clusters
combined.sct@meta.data = combined.sct[[]] %>% 
  mutate(
    integrated_clusters = case_when(
      seurat_clusters %in% c(1,2,5,10,11,13,15,16,19) ~ 'D2',  
      seurat_clusters %in% c(0,3,4,6,7,12,14,18,20) ~ 'D1', 
      seurat_clusters %in% c(8,9,21) ~ 'eSPN', 
      seurat_clusters %in% c(17) ~ 'IC'), 
    
    integrated_subclusters = case_when(
      seurat_clusters %in% c(1,2,5,11,16,19) ~ 'D2-Matrix',  
      seurat_clusters %in% c(10,13,15) ~ 'D2-Striosome',  
      seurat_clusters %in% c(0,3,4,7,12,14,18,20) ~ 'D1-Matrix', 
      seurat_clusters %in% c(6) ~ 'D1-Striosome', 
      seurat_clusters %in% c(8) ~ 'D1/D2H', 
      seurat_clusters %in% c(9,21) ~ 'D1-NUDAP', 
      seurat_clusters %in% c(17) ~ 'IC')
  )

table(is.na(combined.sct$integrated_clusters), combined.sct$Project)
table(is.na(combined.sct$integrated_clusters), combined.sct$seurat_clusters)
table(is.na(combined.sct$integrated_subclusters), combined.sct$seurat_clusters)

## make the plots
pdf(here::here(PLOTDIR, 'Replicates_plots_bySpecies_mainCluster.umap.pdf'), 
    width = 8, height = 5)
DimPlot(combined.sct, reduction = "umap", split.by = "Species", 
        cols = ArchR::paletteDiscrete(combined.sct$integrated_clusters),
        label = T, group.by = 'integrated_clusters', raster = F) + 
  theme(legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 3))) &
  FontSize(main = 7, x.title = 7, y.title = 7)
dev.off()

## make the plots
cols = setNames(RColorBrewer::brewer.pal(7,'Paired'), 
                c('D1-Matrix', 'D1-Striosome', 'D1-NUDAP',  'D1/D2H', 
                  'D2-Matrix', 'D2-Striosome', 'IC'))

combined.sct@meta.data$integrated_subclusters = 
  factor(combined.sct$integrated_subclusters, names(cols))

pdf(here::here(PLOTDIR, 'Replicates_plots_bySpecies_subCluster.umap.pdf'), 
    width = 8, height = 5)
DimPlot(combined.sct, reduction = "umap", split.by = "Species", raster = F,
        cols = cols, label = T, group.by = 'integrated_subclusters') +
  theme(legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 2, override.aes = list(size = 3))) &
  FontSize(main = 7, x.title = 7, y.title = 7)
dev.off()

##################################
## 4) save the integrated object
save_fn = here(DATADIR,'rdas/integrated_mouse_rat_macaque_human_striatum.ARC.rds')
SaveSeuratRds(combined.sct, save_fn)




# ## specify order of integration defined by
# ## https://satijalab.org/seurat/reference/integratedata
# integrate.matrix = matrix(c(-1, -2, # 1: the mice, RNA
#                             -3, 1,  # 2: the rodents, RNA
#                             -4, -5, # 3: the monkeys, RNA
#                             -6, -7, # 4: the humans, RNA
#                             3,  4,  # 5: the primates, outputs of lines 3 & 4
#                             2, -8,  # 6: the rodents, RNA + ATAC
#                             5, -9,   # 7: the primates, RNA + ATAC
#                             6,  7   # 6: the mammals
# ), ncol = 2, byrow = T)

