library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(data.table)
library(future)
library(viridis)
library(cluster)
library(ggh4x)
library(here)
library(tidymodels)
library(broom)
library(dplyr)
library(metap)

# Enable parallelization
plan("multicore", workers = 12)
options(future.globals.maxSize = 100 * 1024^3, future.seed=T)

DATADIR = 'data/tidy_data/'
PLOTDIR = 'figures/multispecies_integration/plots'
TABLDIR = 'figures/multispecies_integration/tables'
dir.create(PLOTDIR, recursive = T)
dir.create(TABLDIR, recursive = T)

###############################################################################
## 1) recalculate the normalized RNA that is supposed to be good for marker genes
save_fn = here(DATADIR,'rdas/integrated_multispecies_striatum_ARC_P15.filtered.rds')
obj = readRDS(save_fn)

obj = obj %>% subset(Assay == 'RNA')
obj = JoinLayers(obj)

###########################################################
## 2) compute the per cell type differential gene table

# for each cell type, find the 1 vs. all marker genes
Idents(obj) = 'integrated_clusters'
clusters = set_names(unique(obj$integrated_clusters)) %>% sort()
markers_df = lapply(clusters, function(x){
   df = obj %>% 
     FindConservedMarkers(ident.1 = x, grouping.var = 'Project', assay = "RNA", 
                          slot = "data", min.cells.group = 3, verbose = T, 
                          logfc.threshold = 0, min.pct = 0, recorrect_umi = F,
                          meta.method = metap::wilkinsonp) %>% 
     rownames_to_column('gene')
  return(df)
  }) %>% bind_rows(.id = 'celltype')


Idents(obj) = 'integrated_subclusters'
subclusters = set_names(unique(obj$integrated_subclusters)) %>% sort()
markers_df2 = lapply(subclusters, function(x){
  df = obj %>% 
    FindConservedMarkers(ident.1 = x, grouping.var = 'Project', assay = "RNA", 
                         slot = "data", min.cells.group = 3, verbose = T, 
                         logfc.threshold = .01, min.pct = .01, recorrect_umi = F,
                         meta.method = metap::wilkinsonp) %>% 
    rownames_to_column('gene')
  return(df)
}) %>% bind_rows(.id = 'celltype')


#######################################################
## 3) adjust the statistics and add consensus metrics

## correct the nominal p-values, calculate the log2FC
markers_df3 = markers_df2 %>% filter(celltype != 'IC') %>% 
  rbind(markers_df) %>% 
  mutate(max_pfdr = p.adjust(max_pval,'fdr')) %>% 
  rowwise() %>% 
  mutate(avg_log2FC = mean(c(`Phan et al._avg_log2FC`, 
                             `Tran et al._avg_log2FC`, 
                             `He, Kleyman et al._avg_log2FC`, 
                             `Stanley et al._avg_log2FC`, 
                             `Saunders et al._avg_log2FC`, 
                             `Savell et al._avg_log2FC`))) %>% 
  arrange(max_pval)

## calculate the consistency of the log2FC sign across datasets
markers_df3$proportion_consistent_fc = markers_df3 %>% 
  dplyr::select(contains('al._avg_log2FC')) %>% 
  apply(2, sign) %>% apply(1, function(x) abs(mean(x)) )

## filter out best conserved marker genes for each cell type
alpha = 0.05
markers_list = markers_df3 %>% 
  dplyr::relocate(c( avg_log2FC, starts_with('max'), proportion_consistent_fc,
                     contains('log2FC'), contains('p_val')), 
                  .after = 'gene') %>% 
  filter(avg_log2FC > 0, max_pfdr < alpha)  %>% 
  dplyr::select(-contains('pct')) %>% 
  split(f = .$celltype)

markers_list = markers_list[c(clusters[-4], subclusters[-7])]
names(markers_list) = make.names(names(markers_list))

########################################
## 3) export the conserved marker genes
here(DATADIR,'rdas/conserved_striatal_neuron_differential_markers2.rds') %>% 
  saveRDS(object = markers_df3)

dir.create(here(DATADIR,'tables'))
here(DATADIR,'tables/conserved_striatal_neuron_marker_genes2.xlsx') %>% 
  writexl::write_xlsx(x = markers_list)

#################################################
# 4) check number of overlapping markers for text
tmp = markers_list %>% lapply(filter, proportion_consistent_fc ==1)
sapply(tmp , nrow)
#     D1           D2         eSPN    D1.Matrix     D1.NUDAP D1.Striosome 
#     43           55           54           37           24           20 
# D1.D2H    D2.Matrix D2.Striosome 
#     57           57           34 

gene_list = tmp %>% lapply('[[', 'gene')

grep('FOXP2', markers_list[['eSPN']]$gene, value = T) ## third top gene
grep('TSHZ1', markers_list[['eSPN']]$gene, value = T) ## third top gene
gene_list[c(1,3)] %>% Reduce(f = 'intersect')  # 0 genes, shared
gene_list[c(2,3)] %>% Reduce(f = 'intersect') # 0 genes shared

## shared markers D1R and D2R onlys
gene_list[c(1,2)] %>% Reduce(f = 'intersect') # 0 gene shared


