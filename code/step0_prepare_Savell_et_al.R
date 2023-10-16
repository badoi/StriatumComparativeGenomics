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

save_fn = here(DATADIR,'Rat_Savell/Savell2020_snRNA_filteredImputed_SeuratObj_N4.h5Seurat')
obj = LoadH5Seurat(save_fn)
head(obj[[]])

obj2 = subset(obj, Stim=='Saline')
counts = obj2@assays$RNA@counts
meta = obj2[[]]

######################################
## 1) create object w/ just human genes

## read in the 1-1 gene orthologs
hg38_to_mm10_genes = read_tsv('/home/bnphan/resources/genomes/GRCh38.p13/ENSEMBL_GRCh38.p13_genes_to_Orthologous_mouse_genes.tsv', show_col_types = FALSE)%>% 
  rename_with(make.names) %>% dplyr::select(c(Mouse.gene.name, Gene.name)) %>% deframe()

counts_hg = counts
rownames(counts_hg) = hg38_to_mm10_genes[rownames(counts)]
counts_hg = counts_hg[!is.na(rownames(counts_hg)) & !duplicated(rownames(counts_hg)), ]
counts_hg[1:5,1:5]

obj3 = CreateSeuratObject(counts_hg, project = "Savell2020", assay = "RNA",
                         meta.data = meta) %>% 
  NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()

table(obj3$CellType)

obj4 = obj3[,grepl('MSN', obj3$CellType)] %>% RunPCA()

save_fn2 = here(DATADIR,'Rat_Savell/Savell2020_rat_striatum_msns.hgGenes.h5Seurat')
SaveH5Seurat(obj4, save_fn2, overwrite = T)





