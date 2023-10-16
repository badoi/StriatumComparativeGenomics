library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(rliger)
library(here)


DATADIR = 'data/tidy_data/'

## loads in sce.nac.tran
load('data/tidy_data/Human_Tran/SCE_NAc-n8_tran-etal.rda')

sce.nac.tran$cellType %>% table()

## subset to just MSNs from Tran et al.
sce.nac.tran.msn = sce.nac.tran[,grepl('MSN', sce.nac.tran$cellType)]
obj_msn = as.Seurat(sce.nac.tran.msn)
save_fn = 'data/tidy_data/Human_Tran/SCE_NAc-n8_tran-etal.msns.h5Seurat'
SaveH5Seurat(obj_msn, save_fn)
