library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(future)
library(here)

# Enable parallelization
plan("multicore", workers = 4)
options(future.globals.maxSize = 100 * 1024^3)

DATADIR = 'data/tidy_data/'

## loads in sce.nac.tran
obj = readRDS(here(DATADIR, 'Mouse_Zeng/WMB-10Xv3-STR.seurat.rds'))

head(obj,2)

table(obj$division)
table(obj$region_of_interest_acronym)
table(obj$class, obj$region_of_interest_acronym)

ind_msn = which(obj$region_of_interest_acronym %in% c('STRd', 'STRv') & 
                obj$class %in% c('06 CNU GABA') & 
                  obj$subclass %in% 
                  c('051 MSN D1 Gaba', '052 MSN D2 Gaba', '050 OT D3 Folh1 Gaba', 
                    '053 MSN D1 Sema5a Gaba', '054 STR-PAL Chst9 Gaba'))

obj_msn = obj[,ind_msn] 

## subset to just MSNs from unaffected individuals
save_fn = here(DATADIR, 'Mouse_Zeng/WMB-10Xv3-STR.msn.seurat.rds')
saveRDS(obj_msn, save_fn)

## convert to human genes
## read in the 1-1 gene orthologs
hg38_to_mm10_genes = read_tsv('/home/bnphan/resources/genomes/GRCh38.p13/ENSEMBL_GRCh38.p13_genes_to_Orthologous_mouse_genes.tsv', show_col_types = FALSE)%>% 
  rename_with(make.names) %>% dplyr::select(c(Mouse.gene.name, Gene.name)) %>% deframe()

counts_hg = LayerData(obj_msn, layer = 'counts')
rownames(counts_hg) = hg38_to_mm10_genes[rownames(counts_hg)]
counts_hg = counts_hg[!is.na(rownames(counts_hg)) & !duplicated(rownames(counts_hg)), ]
counts_hg[1:5,1:5]

obj = CreateSeuratObject(counts_hg, project = "Zeng2023", assay = "RNA",
                         meta.data = obj_msn[[]]) %>% 
  NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  SCTransform(conserve.memory = T, verbose = F)

dir.create(here::here(DATADIR,'Mouse_Zeng'), showWarnings = F)
save_fn2 = here::here(DATADIR,'Mouse_Zeng/Zeng2023_mouse_striatum_msns.hgGenes.rds')
saveRDS(obj, save_fn2)

