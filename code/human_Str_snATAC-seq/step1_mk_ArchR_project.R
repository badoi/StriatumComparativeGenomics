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
addArchRGenome("hg38")
DATADIR='data/raw_data'

################################################
## 0) grab the CATlas annotated cells and samples 
pd = readxl::read_xlsx(
  here('data/raw_data', 'tables', 'Metatable of snATAC-seq experiments.xlsx'), 
  sheet = 'Sheet 2') %>% rename_with(make.names) %>%
  mutate_if(is.factor, as.character)

samples = pd %>% pull(Sample) %>% set_names()

###############################
# 1) create the ArchR project
dir.create(here('data/tidy_data/ArchRProjects'), showWarnings = F)
projectName = here("data/tidy_data/ArchRProjects/human_Str_snATAC")

#create ArchR project out of arrow files
ArrowFiles = here("data/raw_data/arrow/") %>% paste0(samples, '.arrow')
names(ArrowFiles) = samples
all(file.exists(ArrowFiles))
proj <- ArchRProject(ArrowFiles, outputDirectory = projectName, 
                     copyArrows = T)
proj <- saveArchRProject(proj)

#############################
## 2) filter out bad cells ##
proj = filterDoublets( proj, cutEnrich = .5, cutScore = -Inf, filterRatio = 1)

## increase unique fragments cutoff to 10^3, remove cluster of low QC cell 
idxSample <- BiocGenerics::which(proj$nFrags > 10^3)
proj = subsetCells(ArchRProj = proj, cellNames = proj$cellNames[idxSample])
proj <- saveArchRProject(proj)

## fit probability regression for outlier cell detection
df = getCellColData(proj) %>% as.data.frame() %>% mutate(tmp = seq(n())) %>%
  rownames_to_column("CellBarcode")

resid_cutoff = 2
cellsSampleKeep = split(df, df$Sample) %>% 
  map(., ~glm.nb(nFrags ~ TSSEnrichment + PromoterRatio + DoubletEnrichment, 
                 data = .)) %>% 
  map2(.x = ., .y = split(df, f = df$Sample), 
       .f = ~ augment_columns(x = .x, data = .y)) %>% 
  bind_rows() %>% mutate(outlierScore =  abs(.std.resid), 
                         isOutlier = outlierScore > resid_cutoff) %>%
  dplyr::select(-starts_with('\\.')) %>% arrange(tmp) %>%
  filter(!isOutlier) %>% pull(CellBarcode)

## remove cells with low unique fragments
proj = subsetCells(ArchRProj = proj, cellNames = cellsSampleKeep)
proj <- saveArchRProject(proj)

table(proj$Sample)

################################################
## 3) add the sample-level metadata to the cells
pd2 = pd[match(getCellNames(proj) %>% ss('#'), pd$Sample), ] %>% 
  dplyr::select(-Sample)

for(col in names(pd2)){
  proj = addCellColData(proj, data = pd2 %>% pull(col), name = col,
                        cells = getCellNames(proj), force = TRUE)
}

table(proj$Region)

################################################################################
## 4) add the cell-level metadata to the cells from individually annotated cells
proj2 = here('/projects/pfenninggroup/singleCell/BICCN_human_CATlas_snATAC-seq',
             'data/tidy_data/ArchRProjects/BICCN_human_Str_snATAC') %>% 
  loadArchRProject()
proj3 = here('/projects/pfenninggroup/machineLearningForComputationalBiology',
             'snATAC_cross_species_caudate/data/raw_data/hg38',
             'Corces_2020/ArchR_Corces2020_caudate_labeled') %>% 
  loadArchRProject()

pd_cell2 = getCellColData(proj2) %>% data.frame() %>% 
  rownames_to_column('CellID') %>%
  dplyr::select(-c(TSSEnrichment:match, ClustersX300)) %>% 
  rename('PublicationCelltype' = 'celltype')
pd_cell3 = getCellColData(proj3) %>% data.frame() %>% 
  rownames_to_column('CellID') %>%
  dplyr::select(-c(TSSEnrichment:Cluster, is_Dopaminergic., ReadsInPeaks, FRIP, 
                   Cell_Type_Group:Neuron_Type_Group)) %>% 
  rename('PublicationCelltype' = 'FinalClusters')

## gather the MSN annotations from the multi-species integration
pd_msn = here('data/tidy_data/rdas/integrated_multispecies_striatum.meta.rds') %>% 
  readRDS() %>% filter(Assay == 'ATAC', Species == 'Human') %>% 
  rownames_to_column('CellID') %>% dplyr::select(-c(Project:seurat_clusters.full.score))

## combine the prior annotations together 
pd_cell = rbindlist(list(pd_cell2, pd_cell3), fill = T) %>% 
  mutate(Clusters2 = case_when(!grepl('MSN|INT', Clusters2) ~ Clusters2, 
                               TRUE ~ ss(Clusters2, '_')))
table(pd_cell$Sample, pd_cell$Clusters2)

pd_cell = left_join(pd_cell, pd_msn) %>% 
  mutate(integrated_clusters = 
           case_when(integrated_clusters == 'IC' ~ 'D1', 
                     str_detect('MSN', Clusters2) & !is.na(integrated_clusters) ~  integrated_clusters, 
                     str_detect('MSN', Clusters2) & is.na(integrated_clusters) ~ 'Drop', 
                     T ~ Clusters2), 
         integrated_clusters = make.names(integrated_clusters), 
         integrated_subclusters = 
           case_when(integrated_subclusters == 'D1-NUDAP' ~ 'D1/D2H', 
                     integrated_subclusters == 'IC' ~ 'D1-Matrix', 
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
pd_cell2 = pd_cell[ pd_cell$CellID %in% getCellNames(proj) , -c('Sample')]
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

save_dir = here("data/tidy_data/ArchRProjects/human_Str_snATAC_annotated")
proj2 = subsetArchRProject(ArchRProj = proj, cells = cells_keep,
  outputDirectory = save_dir,  dropCells = TRUE, force = TRUE)
