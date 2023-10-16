library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(rliger)
library(future)
library(here)

# Enable parallelization
plan("multicore", workers = 24)
options(future.globals.maxSize = 100 * 1024^3)

DATADIR = 'data/tidy_data/'

## loads in sce.nac.tran
obj_msn = LoadH5Seurat('/projects/pfenninggroup/singleCell/Logan_BU_Striatum_snRNA-seq/data/tidy_data/Seurat_projects/BU_OUD_Striatum_refined_msn_SeuratObj_N22.h5Seurat')

SaveH5Seurat(obj_msn, here(DATADIR, 'Human_Phan/BU_OUD_Striatum_refined_msn_SeuratObj_N22.h5Seurat'))

obj_msn$celltype3 %>% table()
obj_msn$DSM.IV.OUD %>% table()

## subset to just MSNs from unaffected individuals
obj_msn_ctl = obj_msn[,obj_msn$DSM.IV.OUD == 'CTL']
save_fn = here(DATADIR, 'Human_Phan/OUD_Striatum_refined_SeuratObj_N10.msns.h5Seurat')
SaveH5Seurat(obj_msn_ctl, save_fn)
