library(tidyverse)
library(SingleCellExperiment)
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(data.table)
library(future)
library(ArchR)
library(here)

# Enable paralleZuzation
plan("multicore", workers = 24)
options(future.globals.maxSize = 100 * 1024^3)

DATADIR = 'data/tidy_data/'

proj = file.path('/projects/pfenninggroup/singleCell/BICCN_human_CATlas_snATAC-seq',
                 'data/tidy_data/ArchRProjects/BICCN_human_Str_snATAC_MSN') %>% 
  loadArchRProject(showLogo = F)

getCellColData(proj)
table(proj$Clusters2)

mat_se = getMatrixFromProject(proj, useMatrix = 'GeneScoreMatrix')
counts = assays(mat_se)$GeneScoreMatrix
rownames(counts) = rowData(mat_se)$name
counts[1:5,1:5]
meta = colData(mat_se) %>% as.data.frame()

counts[1:5,1:5]
counts = round(counts*10)

obj = CreateAssayObject(counts) %>% 
  CreateSeuratObject(project = "Li2023", assay = "RNA", meta.data = meta) %>% 
  NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()

dir.create(here::here(DATADIR,'Human_Li'), showWarnings = F)
save_fn2 = here::here(DATADIR,'Human_Li/Li2023_human_striatum_msns.rds')
saveRDS(obj, save_fn2)



