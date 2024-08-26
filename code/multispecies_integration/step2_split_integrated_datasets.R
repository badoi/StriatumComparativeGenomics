library(tidyverse)
library(data.table)
library(Seurat)
library(here)

# Enable parallelization
plan("multicore", workers = 12)
options(future.globals.maxSize = 100 * 1024^3, future.seed=TRUE)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

DATADIR = 'data/tidy_data/'

subclusters = c('D1-Matrix', 'D1-Striosome', 'D1-NUDAP',  'D1/D2H', 
                'D2-Matrix', 'D2-Striosome', 'IC')
clusters = c('D1', 'eSPN', 'D2', 'IC')
cols1 = setNames(RColorBrewer::brewer.pal(7,'Paired'), subclusters)
cols2 = setNames(RColorBrewer::brewer.pal(7,'Paired')[c(2,4,6,7)], clusters)

proj1 = c('Gayden et al.','Siletti et al.', 'Tran et al.', 'Phan et al.', 
          'Li et al.', 'Corces et al.', 'Chiou et al.', 'He, Kleyman et al.', 
          'Phillips et al.', 'Savell et al.', 'Zeng et al.', 'Chen et al.',
          'Saunders et al.',  'Stanley et al.',  'Zu et al.')

proj2 = proj1 %>% str_replace('Gayden', 'Gayden Kozel')

# load the metadata of all cells
metadata = here('data/tidy_data/rdas/integrated_multispecies_striatum.meta.rds') %>% 
  readRDS() %>% 
  mutate(integrated_subclusters = factor(integrated_subclusters, levels = subclusters), 
         integrated_clusters = factor(integrated_clusters, levels = clusters),
         Project = as.character(Project), 
         Replicate = ifelse(Project == 'Zeng et al.', ss(CellID, '_', 2), 
                            as.character(Replicate)),
         Replicate = ifelse(Project == 'Siletti et al.', ss(CellID, ':', 1), 
                            Replicate),
         Replicate = ifelse(Project == 'Chiou et al.', ss(CellID, '-', 1), 
                            Replicate),
         Case = case_when(grepl('13114|13151|612|13291|1488|1034', Replicate) & 
                            Project == 'Phan et al.' ~  'Unaffected', T ~ Case), 
         Project = ifelse(is.na(Project), 'Gayden et al.', as.character(Project)),
         Project = factor(Project, levels = proj1))

levels(metadata$Project) = proj2
metadata$Replicate = ifelse(metadata$Project == 'Gayden Kozel et al.', 
                            ss(metadata$CellID, '_', 3), metadata$Replicate)
saveRDS(metadata, here('data/tidy_data/rdas/integrated_multispecies_striatum.meta.rds'))

table(metadata$Project, metadata$Case)
table(metadata$Replicate)
table(metadata$Project)

metadata %>% filter(Project == 'Gayden Kozel et al.') %>% count(Replicate)
metadata %>% filter(Project == 'Gayden Kozel et al.') %>% pull(CellID) %>% head()

###############################################################################
## 1) recalculate the normalized SCT that is supposed to be good for marker genes
save_fn = here(DATADIR,'rdas/integrated_multispecies_striatum_ARC_P15.filtered.rds')
object = readRDS(save_fn)

## make sure not split layers
DefaultAssay(object) = 'RNA'
if (!'counts' %in% Layers(object)) {
  object = JoinLayers(object)
  object %>% saveRDS(save_fn)
}

## update the metadata of this object
object$CellID = Cells(object)
object$integrated_clusters = factor(object$integrated_clusters, levels = clusters)
object$integrated_subclusters = factor(object$integrated_subclusters, levels = subclusters) 
object$Project = as.character(object$Project)

## re-split the replicates to the level of biospecimens/libraries
object$Replicate = metadata$Replicate[match(object$CellID, metadata$CellID)]
object$Case = metadata$Case[match(object$CellID, metadata$CellID)]

## change the levels of the project to update Gayden to Gayden Kozel
object$Project = factor(object$Project, levels = proj1)
levels(object$Project) = proj2
all.equal(Cells(object), rownames(object[[]]))

table(object$Project, object$Case)
table(object$Replicate)

## save the updated object
DefaultAssay(object) = 'RNA'
object = JoinLayers(object)

saveRDS(object, save_fn)

##############################################
## 2) split the integrated object by project
for (project in unique(object$Project)){
  ## subset integrated object by project
  obj_subset = object %>% subset(Project == project)
  
  ## get the species and assay of this project
  species = unique(obj_subset$Species)
  assay = unique(obj_subset$Assay)
  
  save_fn =  here(DATADIR,'rdas',
                  paste0(species,'.',assay,'.',make.names(project), 
                         'integrated_labels.rds'))

  if(!file.exists(save_fn)){
  ## remove metadata columns that have all NA values
  keep = apply(obj_subset[[]], 2, function(x) !all(is.na(x))) %>% which()
  obj_subset@meta.data = obj_subset[[]][,keep]
    
  SaveSeuratRds(obj_subset, save_fn)
}
}

###################################################################
## 3) compute the pseudobulk expression matrix across all datasets


