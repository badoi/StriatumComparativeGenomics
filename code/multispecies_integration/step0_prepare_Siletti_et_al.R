library(tidyverse)
library(Seurat)
library(future)
library(here)

# Enable parallelization
plan("multicore", workers = 4)
options(future.globals.maxSize = 100 * 1024^3)

DATADIR = 'data/tidy_data/'

## loads in sce.nac.tran
obj = readRDS(here(DATADIR, 'Human_Siletti/WHB-10Xv3-Neurons-Basal.nuclei.1.seurat.rds'))

head(obj,2)
table(obj$cluster_alias)
table(obj$description_supercluster)
table(obj$region_of_interest_label)

ind_msn = which(obj$description_supercluster %in% 
                  c('Medium spiny neuron', 'Eccentric medium spiny neuron') & 
                  obj$region_of_interest_label %in% 
                  c('Human CaB', 'Human NAC', 'Human Pu'))

obj_msn = obj[,ind_msn] 
# obj_msn = obj_msn %>% SCTransform(conserve.memory = TRUE, verbose = T)

## subset to just MSNs from unaffected individuals
save_fn = here(DATADIR, 'Human_Siletti/WHB-10Xv3-Neurons-Basal.nuclei.msn.seurat.rds')
saveRDS(obj_msn, save_fn)
