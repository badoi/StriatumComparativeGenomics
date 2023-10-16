library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(data.table)
library(future)
library(viridis)
library(cluster)
library(ggh4x)
library(here)
library(tidymodels)
library(broom)
library(dplyr)
library(metap)

# Enable parallelization
plan("multicore", workers = 12)
options(future.globals.maxSize = 100 * 1024^3, future.seed=TRUE)

DATADIR = 'data/tidy_data/'

###############################################################################
## 1) recalculate the normalized SCT that is supposed to be good for marker genes
save_fn = here(DATADIR,'rdas/integrated_mouse_rat_macaque_human_striatum.ARC.h5Seurat')
obj = LoadH5Seurat(save_fn)

for (project in unique(obj$Project)){
  ## subset integrated object by project
  obj_subset = obj %>% subset(Project == project)
  
  ## get the species and assay of this project
  species = unique(obj_subset$Species)
  assay = ifelse(project %in% c('Corces et al.', 'Li et al.'), 'ATAC', 'RNA')
  
  save_fn =  here(DATADIR,'rdas',paste0(species,'.',assay,'.',make.names(project), 'integrated_labels.h5Seurat'))
  if(!file.exists(save_fn)){
   
    ## remove metadata columns that have all NA values
    keep = apply(obj_subset[[]], 2, function(x) !all(is.na(x))) %>% which()
    obj_subset@meta.data = obj_subset[[]][,keep]
      
    SaveH5Seurat(obj_subset, save_fn)
  }
}