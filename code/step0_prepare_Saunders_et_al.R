library(DropSeq.util)
library(SingleCellExperiment)
library(tidyverse)
library(scater)
library(BiocParallel)
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(data.table)
library(future)
library(here)

# Enable parallelization
plan("multicore", workers = 24)
options(future.globals.maxSize = 100 * 1024^3, stringsAsFactors = F)

DATADIR = 'data/tidy_data/'
DATA2= '/projects/pfenninggroup/singleCell/Saunder_et_al_DropViz'

# DropViz metadata
annotation = readRDS(here(DATA2, 'data/annotation.BrainCellAtlas_Saunders_version_2018.04.01.RDS')) %>% filter(tissue == 'STR')

# add cluster names
clusters = read.csv(here(DATA2, 'data/clusters.csv'))
clusters = clusters[clusters$Region %in% c('Striatum'),]

####################################
# 1) paths to the striatum data
counts.path <- here(DATA2, "data/F_GRCm38.81.P60Striatum.raw.dge.txt.gz")
meta.path <- here(DATA2,"data/F_GRCm38.81.P60Striatum.cell_cluster_outcomes.RDS")

# load in the data
counts <- loadSparseDge(counts.path) 
meta = readRDS(meta.path)
indSTR = match(paste0('STR_',meta$subcluster), annotation$tissue_subcluster)
meta = cbind(meta, annotation[indSTR, names(annotation) != 'subcluster'])

table(meta$reason)
obj = CreateSeuratObject(counts, project = "Saunders_2018", assay = "RNA",
                         meta.data = meta) 
table(obj$reason)

obj = obj[,is.na(obj$reason)]%>% 
  NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()

save_fn = here(DATADIR,'Mouse_Saunders/Saunders2018_mouse_striatum_all.h5Seurat')
SaveH5Seurat(obj, save_fn)


## create object w/ just human genes

##########################################
## 2) convert to human 1-1 gene orthologs
hg38_to_mm10_genes = read_tsv('/home/bnphan/resources/genomes/GRCh38.p13/ENSEMBL_GRCh38.p13_genes_to_Orthologous_mouse_genes.tsv', show_col_types = FALSE)%>% 
  rename_with(make.names) %>% dplyr::select(c(Mouse.gene.name, Gene.name)) %>% deframe()

counts_hg = obj@assays$RNA@counts
meta = obj[[]]
rownames(counts_hg) = hg38_to_mm10_genes[rownames(counts)]
counts_hg = counts_hg[!is.na(rownames(counts_hg)) & !duplicated(rownames(counts_hg)), ]
counts_hg[1:5,1:5]

obj2 = CreateSeuratObject(counts_hg, project = "Saunders2018", assay = "RNA",
                         meta.data = meta) %>% 
  NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()

meta %>% filter(class == 'NEURON') %>% count(common_name)
obj3 = obj2[, obj2$class == 'NEURON' & grepl('SPN', obj2$common_name)] %>% 
  FindVariableFeatures() %>% ScaleData() %>% RunPCA()
table(obj3$common_name)
obj3$common_name = str_replace_all(obj3$common_name, 'eccentric ', 'e')

save_fn2 = here(DATADIR,'Mouse_Saunders/Saunders2018_mouse_striatum_msns.hgGenes.h5Seurat')
SaveH5Seurat(obj3, save_fn2, overwrite = T)





