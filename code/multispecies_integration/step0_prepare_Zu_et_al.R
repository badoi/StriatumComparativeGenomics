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

proj = file.path('/projects/pfenninggroup/singleCell/BICCN_mouse_CATlas_snATAC-seq',
                 'data/tidy_data/ArchRProjects/BICCN_mouse_Str_snATAC2_MSN') %>% 
  loadArchRProject(showLogo = F)

getCellColData(proj)
table(proj$Clusters2)

mat_se = getMatrixFromProject(proj, useMatrix = 'GeneScoreMatrix')
counts = assays(mat_se)$GeneScoreMatrix
rownames(counts) = rowData(mat_se)$name
counts[1:5,1:5]
meta = colData(mat_se) %>% as.data.frame()

######################################
## 1) create object w/ just human genes

## read in the 1-1 gene orthologs
hg38_to_mm10_genes = read_tsv('/home/bnphan/resources/genomes/GRCh38.p13/ENSEMBL_GRCh38.p13_genes_to_Orthologous_mouse_genes.tsv', show_col_types = FALSE)%>% 
  rename_with(make.names) %>% dplyr::select(c(Mouse.gene.name, Gene.name)) %>% deframe()

counts_hg = counts
rownames(counts_hg) = hg38_to_mm10_genes[rownames(counts)]
counts_hg = counts_hg[!is.na(rownames(counts_hg)) & !duplicated(rownames(counts_hg)), ]
counts_hg[1:5,1:5]
counts_hg = round(counts_hg*10)

obj = CreateSeuratObject(counts_hg, project = "Zu2023", assay = "RNA",
                         meta.data = meta) %>% 
  NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()

dir.create(here::here(DATADIR,'Mouse_Zu'), showWarnings = F)
save_fn2 = here::here(DATADIR,'Mouse_Zu/Zu2023_mouse_striatum_msns.hgGenes.rds')
saveRDS(obj, save_fn2)



