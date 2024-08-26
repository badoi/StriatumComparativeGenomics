## conda activate r4
## packages for data table processing 
library(here)
library(tidyverse)
library(parallel)

## main Seurat package snRNA-seq pacakges
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(DropletQC)
library(future)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

###################################################################
# 0) Seurat uses the future package for parallelization
## set to be parallel over 28 cores with 500Gb
plan("multicore", workers = 28)
options(future.rng.onMisuse = 'ignore')
options(future.globals.maxSize = 150 * 1024^3)

###########################################################################
# 1) load in indvidual snRNA-seq objects and create merged Seurat projects
save_fn = list.files(here('data/raw_data/Seurat_objects'),
                     pattern = 'STARsolo_SoupX_rawCounts', full.names = T) %>% 
  str_subset('25', negate = T)
names(save_fn) = basename(save_fn) %>% ss('STARsolo_SoupX_rawCounts_|.rds', 2)
num_samples = length(save_fn)

objList = mclapply(save_fn, readRDS, mc.cores = 4)
obj = merge(objList[[1]], objList[-1], add.cell.ids = names(objList)) %>% JoinLayers()
obj[['RNA']] = split(obj[['RNA']], obj$orig.ident)
obj = obj %>% ScaleData(verbose = F) %>% RunPCA(verbose = F)

########################################################
# 2) use Seurat reciprocal PCA to join samples together
## find integrating features
obj = IntegrateLayers(object = obj, method = RPCAIntegration, 
  orig.reduction = "pca", new.reduction = "integrated.rpca", k.anchor = 25,
  normalization.method = "LogNormalize", dims = 1:30, verbose = F) 

obj = obj %>% RunUMAP(reduction = "integrated.rpca", dims = 1:30) %>% 
  FindNeighbors(reduction = "integrated.rpca", dims = 1:30, verbose = T) %>% 
  FindClusters(resolution = 1, algorithm = 2, verbose = T)

obj = JoinLayers(obj)

## clean up RAM
rm(objList); gc()

#####################################
# 3) add in patient/sample metadata
pheno = here('data/tidy_data/tables/Gayden_NAc_multiome_Sample_sheet.xlsx') %>% 
  readxl::read_xlsx() %>% rename_with(make.names) %>%
  rename_with(~ gsub("(\\.){2,}", '\\.', .x)) %>%
  rename_with(~ gsub('\\.$', '', .x)) %>% column_to_rownames('GEX')

## take a look at the phenotype table
head(pheno)

## look at the per-cell meta data
head(obj@meta.data)

## append pt. phenotypes to single cell meta data 
obj@meta.data = cbind(obj@meta.data, pheno[obj[[]][,'orig.ident'],])
head(obj@meta.data)

###############################################
# 5) filter the lower quality cells

## Run DropletQC's identify_empty_drops function
nf.umi <- obj[[]] %>% mutate(nf=dropletQC.nucFrac, umi=nCount_RNA) %>% 
  relocate(nf, umi, .before = everything())

## estimate cell, empty droplet, and damaged cell w/ DropletQC algorithms
DropletQC.ed <- nf.umi %>% 
  identify_empty_drops(nf_rescue = 0.50, umi_rescue = 1000) %>%
  relocate(cell_status, seurat_clusters, .after = umi) %>%
  dplyr::select(nf:seurat_clusters) %>%
  identify_damaged_cells(nf_sep = 0.15, umi_sep_perc = 50, verbose = F)
DropletQC.ed = DropletQC.ed$df

## add Droplet QC empty droplet estimation to metadata
obj$dropletQC.keep = DropletQC.ed[colnames(obj), 'cell_status']

## look at which clusters should be kept by miQC per-cell fraction
obj@meta.data %>% group_by(seurat_clusters) %>%
  summarise(num = n(), prop = sum(miQC.keep == 'keep') / n() ) %>% 
  arrange(prop)

## look at which clusters should be kept by doublet SCDS per-cell fraction
obj@meta.data %>% group_by(seurat_clusters) %>%
  summarise(num = n(), prop = sum(scds.keep == 'keep') / n() ) %>% 
  arrange(prop)

## look at which clusters should be kept by both metrics
(t1 = obj@meta.data %>% group_by(seurat_clusters) %>%
    summarise(num = n(), 
              numKeep = sum(scds.keep == 'keep' & miQC.keep == 'keep' & 
                              dropletQC.keep == 'cell'), 
              prop = sum(numKeep) / n() ) %>% arrange(prop))

## keep cells in the clusters that have more than 10% of OK cells
good_clusters <- t1 %>% filter(prop > 0.10) %>% pull(seurat_clusters)

## export unfiltered per-cell QC metric table
obj@meta.data = obj[[]] %>% relocate(dropletQC.keep, 
                                                   .after = 'dropletQC.nucFrac')
save_qcTale_fn = here('data/tidy_data/tables', 
                      paste0("Gayden_NAc_unfiltered_QC_table_N",num_samples,'.txt.gz'))
write_tsv(obj@meta.data, save_qcTale_fn)

obj$qc_keep = obj$ miQC.keep == "keep" & 
  obj$scds.keep == "keep" & obj$dropletQC.keep == 'cell' &
  obj$seurat_clusters %in% good_clusters

table(obj$qc_keep, obj$orig.ident)


##################################################################
# 6) save normalized, UMAP embedded, object for downstream analyses
dir.create(here('data/tidy_data/Human_Gayden'), showWarnings = F)
save_proj_fn = here('data/tidy_data/Human_Gayden', 
                    paste0("Gayden_NAc_unfiltered_SCT_SeuratObj_N",num_samples,'.rds'))
saveRDS(obj, save_proj_fn)


## subset cells to those not predicted low QC or doublet
obj_filtered = subset(obj, subset = miQC.keep == "keep" &
                        scds.keep == "keep" & dropletQC.keep == 'cell' &
                        seurat_clusters %in% good_clusters)

## recompute PCA and UMAP embedding post-filtering
obj_filtered = obj_filtered %>% RunPCA(verbose = F) %>%
  RunUMAP(reduction = "integrated.rpca", dims = 1:30) %>%
  FindNeighbors(dims = 1:30, verbose = T) %>%
  FindClusters(resolution = 2, algorithm = 2, verbose = T)

dir.create(here('data/tidy_data/Human_Gayden'), showWarnings = F)
save_proj_fn = here('data/tidy_data/Human_Gayden', 
                    paste0("Gayden_NAc_filtered_SCT_SeuratObj_N",num_samples,'.rds'))
saveRDS(obj_filtered, save_proj_fn)

