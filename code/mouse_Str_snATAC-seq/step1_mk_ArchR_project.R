library(ArchR)
library(parallel)
library(tidyverse)
library(here)
library(MASS)
library(broom)
library(tidymodels)

## add general functions at the top of the R scripts
ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

options(repr.plot.width=15, repr.plot.height=8.5)
addArchRThreads(threads = 12) 
addArchRGenome("mm10")
DATADIR='data/raw_data'

################################################
## 0) grab the CATlas annotated cells and samples 
pd = read_tsv(
  here(DATADIR, 'tables', 'SI_Table_1_Sample_and_dissection_summary.txt')) %>% 
  rename_with(make.names) %>% rename('Sample' = 'sample') %>% 
  mutate(match = Sample, Sample = paste(SubRegion, Sample, sep = '.')) %>% 
  filter(SubRegion %in% c('ACB', 'CP'))

samples = pd %>% pull(Sample) %>% set_names()

###############################
# 1) create the ArchR project
dir.create(here('data/tidy_data/ArchRProjects'), showWarnings = F)
projectName = here("data/tidy_data/ArchRProjects/BICCN_mouse_Str_snATAC2")

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

## increase unique fragments cutoff to 10^3.5, remove cluster of low QC cell 
idxSample <- BiocGenerics::which(proj$nFrags > 10^3.5)
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

################################################
## 3) add the sample-level metadata to the cells
pd2 = pd[match(getCellNames(proj) %>% ss('#'), pd$Sample), ] %>% 
  dplyr::select(-Sample)
for(col in names(pd2)){
  proj = addCellColData(proj, data = pd2 %>% pull(col), name = col,
                        cells = getCellNames(proj), force = TRUE)
}

################################################
## 4) add the cell-level metadata to the cells
pd_cell = read_tsv(
  here(DATADIR, 'tables', '2023-03-05453B-s3/SI_Tables',
       'SI Table 2 Metadata table for all the 2.3 million nuclei in the snATAC-seq data.txt'))  %>% 
  rename('match' = 'Sample') %>% inner_join(x= pd) %>% 
  mutate(CellID = paste0(Sample, '#',Barcode))

table(pd_cell$CellID %in% getCellNames(proj))
table(pd_cell$CellID2 %in% getCellNames(proj))

pd_cell2 = pd_cell[ pd_cell$CellID %in% getCellNames(proj) ,]
pd_cell2 = pd_cell2[, !names(pd_cell2) %in% (getCellColData(proj) %>% names())] %>% 
  dplyr::select(-c(`# of Fragments`, TSSe, Barcode))

table(pd_cell2$Subclass)
table(pd_cell2$NeuronTransmitter)

table(proj$Sample, proj$MajorRegion)

for(col in names(pd_cell2)[-1]){
  proj = addCellColData(proj, data = pd_cell2 %>% pull(col), name = col,
    cells = pd_cell2$CellID, force = TRUE)
}
table(proj$Sample, proj$Subclass)

proj <- saveArchRProject(proj)
