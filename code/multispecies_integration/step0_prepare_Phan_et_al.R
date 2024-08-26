library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(future)
library(here)

# Enable parallelization
plan("multicore", workers = 24)
options(future.globals.maxSize = 100 * 1024^3)

DATADIR = 'data/tidy_data/'

## loads in sce.nac.tran
obj_msn = LoadH5Seurat(here(DATADIR, 'Human_Phan/BU_OUD_Striatum_refined_msn_SeuratObj_N22.h5Seurat'))

obj_msn$celltype3 %>% table()
obj_msn$DSM.IV.OUD %>% table()

## subset to just MSNs from unaffected individuals
save_fn = here(DATADIR, 'Human_Phan/OUD_Striatum_refined_SeuratObj_N22.msns.rds')
saveRDS(obj_msn, save_fn)
