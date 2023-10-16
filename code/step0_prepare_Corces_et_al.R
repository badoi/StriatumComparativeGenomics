library(ArchR)
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(tidyverse)
library(rliger)
library(future)
library(here)

# Enable parallelization
plan("multicore", workers = 24)
options(future.globals.maxSize = 100 * 1024^3)

DATADIR = 'data/tidy_data/'

##################################
## 1 ) loads in the ArchR project
proj = file.path('/projects/pfenninggroup/singleCell/Corces2020_human_snATAC-seq',
                 'data/ArchRProjects/ArchR_Corces2020_human_brain_scATAC-seq') %>% 
  loadArchRProject(showLogo = F)

getCellColData(proj)

mat_se = getMatrixFromProject(proj, useMatrix = 'GeneScoreMatrix')
mat_se = mat_se[,which(mat_se$Region=='CAUD')]
mat_se = mat_se[,which(mat_se$ClusterName %in% c('InhibitoryNeurons'))]

table(meta$ClusterName)

#################################################
## 3 ) extract just the geneScoreMatrix for these
counts = assays(mat_se)$GeneScoreMatrix
counts = round(counts * 10)
rownames(counts) = rowData(mat_se)$name
meta = colData(mat_se) %>% as.data.frame()
counts[1:5,1:5]

obj = CreateSeuratObject(counts, project = "Corces2020", assay = "RNA",
                         meta.data = meta) %>% 
  NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()

dir.create(here(DATADIR,'Human_Corces'), showWarnings = F)
save_fn2 = here(DATADIR,'Human_Corces/Corces2020_human_striatum_msns.hgGenes.h5Seurat')
SaveH5Seurat(obj, save_fn2, overwrite = T)



