library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(data.table)
library(future)
library(here)

# Enable parallelization
plan("multicore", workers = 24)
options(future.globals.maxSize = 100 * 1024^3)

DATADIR = 'data/tidy_data/'

meta = fread(here(DATADIR,'Mouse_Stanley/metadata_final.csv')) %>% 
  column_to_rownames('V1')

counts = fread(here(DATADIR,'Mouse_Stanley/counts.csv')) %>% 
  column_to_rownames('V1') %>% as.matrix()
counts = counts[,rownames(meta)]
counts[1:5,1:5]


obj = CreateSeuratObject(counts, project = "Stanley2020", assay = "RNA",
                         meta.data = meta) %>% 
  NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()

table(obj$discrete.type)
table(obj$continuous.subtype)
save_fn = here(DATADIR,'Mouse_Stanley/Stanley2020_mouse_striatum_msns.rds')
saveRDS(obj, save_fn)


## create object w/ just human genes


## read in the 1-1 gene orthologs
hg38_to_mm10_genes = read_tsv('/home/bnphan/resources/genomes/GRCh38.p13/ENSEMBL_GRCh38.p13_genes_to_Orthologous_mouse_genes.tsv', show_col_types = FALSE)%>% 
  rename_with(make.names) %>% dplyr::select(c(Mouse.gene.name, Gene.name)) %>% deframe()

counts_hg = counts
rownames(counts_hg) = hg38_to_mm10_genes[rownames(counts)]
counts_hg = counts_hg[!is.na(rownames(counts_hg)) & !duplicated(rownames(counts_hg)), ]
counts_hg[1:5,1:5]

obj2 = CreateSeuratObject(counts_hg, project = "Stanley2020", assay = "RNA",
                         meta.data = meta) %>% 
  NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()

table(obj2$discrete.type)
table(obj2$continuous.subtype)
save_fn2 = here(DATADIR,'Mouse_Stanley/Stanley2020_mouse_striatum_msns.hgGenes.rds')
saveRDS(obj2, save_fn2)
