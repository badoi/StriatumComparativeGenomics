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
obj = readRDS(here(DATADIR, 'Macaque_Chiou/Chiou2023_macaque_sciRNA.basal.ganglia.rds'))
head(obj)
head(Features(obj))

table(obj$cell_class)
table(obj$region_name)

ind_msn = which(obj$cell_type %in% c('medium spiny neuron'))
obj_msn = obj[,ind_msn] %>% SCTransform()

## subset to just MSNs from unaffected individuals
save_fn = here(DATADIR, 'Macaque_Chiou/Chiou2023_macaque_sciRNA_basal_ganglia.msn.rds')
saveRDS(obj_msn, save_fn)
