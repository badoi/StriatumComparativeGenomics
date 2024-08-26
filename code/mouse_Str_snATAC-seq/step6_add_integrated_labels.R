library(ArchR)
library(parallel)
library(tidyverse)
library(here)
library(MASS)
library(broom)
library(tidymodels)
library(data.table)

## add general functions at the top of the R scripts
ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

options(repr.plot.width=15, repr.plot.height=8.5)
addArchRThreads(threads = 12) 
addArchRGenome("mm10")
DATADIR='data/raw_data'

###############################
# 1) create the ArchR project
dir.create(here('data/tidy_data/ArchRProjects'), showWarnings = F)
projectName = here("data/tidy_data/ArchRProjects/mouse_Str_snATAC")

## 4) add the cell-level metadata to the cells from individually annotated cells
proj = loadArchRProject(projectName)
table(proj$Clusters2, proj$ClustersX300)

###########################################
# relabel the excitatory and "C4" clusters as low quality cells
remapClust <- c( 'C1' = 'INT', 'C2' = 'INT', 'C3' = 'INT', 'C4' = 'Drop',
  'C5'  = 'INT', 'C6' = 'Drop', 'C7' = 'Drop', 'C8' = 'Drop', 'C9'  = 'Drop', 
  'C10'= 'MSN', 'C11' = 'MSN', 'C12' = 'MSN', 'C13' = 'MSN', 'C14' = 'MSN',
  'C15'  = 'MSN', 'C16' = 'MSN', 'C17' = 'MSN', 'C18' = 'MSN', 'C19'  = 'MSN', 
  'C20'= 'MSN', 'C21' = 'Microglia', 'C22' = 'Oligo', 'C23' = 'OPC',
  'C24' = 'Astro',  'C25' = 'Astro')
proj$Clusters2 <- mapLabels(proj$ClustersX300, newLabels = remapClust, 
                            oldLabels = names(remapClust))
table(proj$Clusters2)

################################################################################
## 2) add the cell-level metadata to the cells from individually annotated cells

## gather the MSN annotations from the multi-species integration
pd_msn = here('data/tidy_data/rdas/integrated_multispecies_striatum.meta.rds') %>% 
  readRDS() %>% filter(Assay == 'ATAC', Species == 'Mouse') %>% 
  rownames_to_column('CellID') %>% dplyr::select(-c(Project:seurat_clusters.full.score))
table(pd_msn$integrated_clusters)
table(pd_msn$integrated_subclusters)

## combine the prior annotations together 
pd_cell = getCellColData(proj) %>% as.data.frame() %>% 
  rownames_to_column('CellID') %>% dplyr::select(CellID, Sample, Clusters2) 

pd_cell = left_join(pd_cell, pd_msn) %>% 
  mutate(integrated_clusters = 
           case_when(integrated_clusters == 'IC' ~ 'Drop', 
                     str_detect('MSN', Clusters2) & !is.na(integrated_clusters) ~  integrated_clusters, 
                     str_detect('MSN', Clusters2) & is.na(integrated_clusters) ~ 'Drop', 
                     T ~ Clusters2), 
         integrated_clusters = make.names(integrated_clusters), 
         integrated_subclusters = 
           case_when(integrated_subclusters == 'IC' ~ 'Drop', 
                     str_detect('MSN', Clusters2) & !is.na(integrated_subclusters) ~ integrated_subclusters,
                     str_detect('MSN', Clusters2) & is.na(integrated_subclusters) ~ 'Drop', 
                     T ~ Clusters2), 
         integrated_subclusters = case_when(integrated_subclusters == 'MSN' ~ 'Drop', 
                                            T ~ integrated_subclusters), 
         integrated_subclusters = make.names(integrated_subclusters))

table(pd_cell$Sample, pd_cell$integrated_clusters)
table(pd_cell$Sample, pd_cell$integrated_subclusters)

## add the cell-level metadata to the big striatum snATAC project
table(pd_cell$CellID %in% getCellNames(proj))
pd_cell2 = pd_cell[ pd_cell$CellID %in% getCellNames(proj) , -c(2,3)]
table(pd_cell2$integrated_subclusters)

for(col in names(pd_cell2)[-1]){
  proj = addCellColData(proj, data = pd_cell2 %>% pull(col), name = col,
    cells = pd_cell2$CellID, force = TRUE)
}

table(proj$Sample, proj$integrated_clusters)
table(proj$Sample, proj$integrated_subclusters)
table(proj$Sample, !is.na(proj$integrated_subclusters))

proj <- saveArchRProject(proj)

## save a subset without dropped cells or NA cells 
cells_keep = getCellColData(proj) %>% data.frame() %>% rownames_to_column('CellID') %>% 
  filter(!is.na(integrated_clusters), integrated_clusters != 'Drop', 
         !is.na(integrated_subclusters), integrated_subclusters != 'Drop') %>% 
  pull(CellID)

save_dir = here("data/tidy_data/ArchRProjects/mouse_Str_snATAC_annotated")
proj2 = subsetArchRProject(ArchRProj = proj, cells = cells_keep,
  outputDirectory = save_dir,  dropCells = TRUE, force = TRUE)
