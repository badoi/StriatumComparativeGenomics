## conda activate r4
## packages for data table processing 
library(here)
library(tidyverse)

## main Seurat package snRNA-seq pacakges
library(Seurat)
library(SeuratDisk)
library(future)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

DATADIR=file.path('/projects/pfenninggroup/singleCell/Logan_BU_Striatum_snRNA-seq', 
                  'data/tidy_data/HeKleyman2021_macaque_striatum_data_processing')

#######################################################
# 0) Seurat uses the future package for parallelization
## set to be parallel over 28 cores
plan("multicore", workers = 28)
options(future.globals.maxSize = 100 * 1024^3)
options(future.rng.onMisuse = 'ignore')

##################################################
# 1) load in cell type labels for label transfer
## read in Logan snRNA dataset to label transfer
save_merged_fn = here('data/tidy_data/Human_Gayden',
                      "Gayden_NAc_filtered_SCT_SeuratObj_N4.rds")
obj_merged = save_merged_fn %>% readRDS()
obj_merged$celltype1 = obj_merged$celltype2 = NULL

## read in labeled monkey MSN subtype dataset, He, Kleyman et al. 2021
## subject to same region in monkey, recompute PCA space
monkey_all = here(DATADIR, 'rdas', 
                  'GSE167920_Results_full_nuclei_processed_final.rds') %>% 
  readRDS() %>% subset(subset = region_name %in% c('nacc')) %>% 
  RunPCA(verbose = F) 
DefaultAssay(monkey_all) = 'integrated'

head(monkey_all[[]])

## dorsal striatum subtypes
monkey_msn = here(DATADIR, 'rdas', 
                  'GSE167920_Results_MSNs_processed_final.v5.rds') %>% 
  readRDS() %>% subset(subset = region_name %in% c( "nacc")) %>% 
  RunPCA(verbose = F) 

DefaultAssay(monkey_msn) = 'RNA'

head(monkey_msn[[]])

###############################################################
# 2) transfer cell type labels from Macaque all nuclei dataset
## compute anchor cells b/t dataset and Macaque  dataset
anchors.mac_all <- FindTransferAnchors(
  reference = monkey_all, query = obj_merged, reduction = 'rpca',
  recompute.residuals = F,
  query.assay = 'RNA', reference.assay = 'integrated', k.anchor = 25,
  normalization.method = 'SCT', dims = 1:30, reference.reduction = "pca")

## predict cell type w/ macaque 'cell_type_2' column
predictions.mac_all <- TransferData(
  anchorset = anchors.mac_all, refdata = monkey_all$cell_type_2, 
  dims = 1:30, store.weights = F)

## append the predicted cell type to 
obj_merged <- AddMetaData(obj_merged, metadata = predictions.mac_all$predicted.id,
                          col.name = 'celltype1')

## see which Seurat clusters corresponds to macaque cell types
(tbl1 = with(obj_merged[[]], table(celltype1, seurat_clusters)))
cluster_to_macAll = colnames(tbl1)
names(cluster_to_macAll) = rownames(tbl1)[apply(tbl1, 2, which.max)]
cluster_to_macAll


######################################################################
# 3) transfer cell type labels from Macaque MSN sub-type nuclei dataset
## split the data to just predicted MSN subtypes
## do this with pipes in one go
obj_msn = obj_merged %>% 
  subset(subset = celltype1 == 'MSNs' | 
           seurat_clusters %in% cluster_to_macAll[names(cluster_to_macAll) =='MSNs'])

obj_msn[['RNA']] = split(obj_msn[['RNA']], obj_msn$orig.ident)
obj_msn = obj_msn %>% FindVariableFeatures(verbose = F) %>% 
  ScaleData(verbose = F) %>% RunPCA(verbose = F)

obj_msn = IntegrateLayers(object = obj_msn, method = RPCAIntegration, k.anchor = 25,
                      orig.reduction = "pca", new.reduction = "integrated.rpca", 
                      normalization.method = "LogNormalize", dims = 1:30, verbose = F) 

obj_msn = obj_msn %>% RunUMAP(reduction = "integrated.rpca", dims = 1:30) %>% 
  FindNeighbors(reduction = "integrated.rpca", dims = 1:30, verbose = F) %>% 
  FindClusters(resolution = 1, algorithm = 2, verbose = F)

## reciprocal PCA project b/t monkey/human MSNs
anchors.mac_msn <- FindTransferAnchors(
  reference = monkey_msn, query = obj_msn, reduction = 'rpca',
  query.assay = 'RNA', reference.assay = 'RNA', k.anchor = 25,
  dims = 1:30, reference.reduction = "pca")

## predict MSNs at the subtype level
predictions.mac_msn <- TransferData(anchorset = anchors.mac_msn, 
                                    refdata = monkey_msn$MSN_type, dims = 1:30)

## transfer MSN subtypes labels back to MSN and full dataset
obj_msn <- AddMetaData(obj_msn, metadata = predictions.mac_msn$predicted.id,
                       col.name = 'celltype2')
obj_merged$celltype2 = obj_msn[['celltype2']][colnames(obj_merged),]
obj_merged$celltype2 = with(obj_merged[[]], ifelse(is.na(celltype2 ),celltype1, celltype2 ))

## see which cluster clusters correponds to MSNs
(tbl2 = with(obj_msn[[]], table(celltype2, seurat_clusters)))
cluster_to_macMSN = colnames(tbl2)
names(cluster_to_macMSN) = rownames(tbl2)[apply(tbl2, 2, which.max)]
table(obj_msn$celltype2)
cluster_to_macMSN

############################################
# 4) save the subset of cells that are MSNs
save_subset_msn = here('data/tidy_data/Human_Gayden',
                       "Gayden_NAc_subsetMSN_SCT_SeuratObj_N4.rds")
saveRDS(obj_msn, file = save_subset_msn)
saveRDS(obj_merged, file = save_merged_fn)

