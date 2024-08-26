library(ArchR)
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(here)
set.seed(1)
ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(repr.plot.width=15, repr.plot.height=8.5)

addArchRThreads(threads = 8) 
addArchRGenome("mm10")

proj <- loadArchRProject(here('data/tidy_data/ArchRProjects', 
                              'BICCN_mouse_Str_snATAC2_MSN'))

varFeats = c(110, 220, 330)

for (varFeat in varFeats){
    iterLSIName = paste0("IterativeLSI",varFeat)
    proj <- addIterativeLSI(
      ArchRProj = proj, 
      clusterParams = list(
        resolution = 0.2, 
        sampleCells = 50000,
        n.start = 10
      ),
      iterations = 4, # increase this if noticing subtle batch effects
      scaleTo = 20000, # median unique fragment per cell
      selectionMethod = 'var',
      sampleCellsPre = 20000,
      varFeatures = varFeat*1000,
      saveIterations = FALSE,
      useMatrix = "TileMatrix",
      depthCol = "nFrags",
      name = iterLSIName,
      force = TRUE
    )

    proj <- saveArchRProject(proj)

    #batchcorrection
    HarmonyName = paste0("HarmonyX",varFeat)
    proj <- addHarmony(proj, reducedDims = iterLSIName,
                       name = HarmonyName, groupBy = 
                         c("Sample", "ProtocalVersion" ,"RegionName", 'SubRegion'), 
                       lambda = rep(1, 4),
                       force = TRUE)
    #clustering
    ClustersName = paste0("ClustersX",varFeat)
    proj <- addClusters(proj, reducedDims = HarmonyName, name = ClustersName, 
                        resolution = .5, maxClusters = 35, force = TRUE)

    #Add UMAP
    UMAPName = paste0("UMAPH",varFeat)
    proj <- addUMAP(
      ArchRProj = proj, reducedDims = HarmonyName, name = UMAPName, 
      nNeighbors = 30, minDist = 0.5, metric = "cosine", force = TRUE)
    proj <- addImputeWeights(proj, reducedDims = HarmonyName)
    proj = saveArchRProject(proj)
}

## add the single cell type clustering from Gayden et al. Striatum mapping
annotated_fn = here::here('/projects/pfenninggroup',
                          'machineLearningForComputationalBiology',
                          'StriatumComparativeGenomics/data/tidy_data/rdas',
                          'Mouse.ATAC.Zu.et.al.integrated_labels.rds')
obj = LoadSeuratRds(annotated_fn)
pd_cell = obj[[]] %>% rownames_to_column('CellID') %>%
  dplyr::select(CellID, integrated_clusters, integrated_subclusters) %>% 
  mutate_if(is.factor, as.character)

for(col in names(pd_cell)[-1]){
  proj = addCellColData(proj, data = pd_cell %>% pull(col), name = col,
                        cells = pd_cell$CellID, force = TRUE)
}
table(proj$Sample, proj$integrated_clusters)
table(proj$Sample, proj$integrated_subclusters)

proj <- saveArchRProject(proj)

