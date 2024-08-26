library(BSgenome.Hsapiens.UCSC.hg38)
library(tidyverse)
library(ArchR)
library(here)

set.seed(1)
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(repr.plot.width=15, repr.plot.height=8.5)

addArchRThreads(threads = 12) 
addArchRGenome("hg38")

############################################
## 1) load the ArchR project and call peaks
save_dir = here("data/tidy_data/ArchRProjects/human_Str_snATAC_annotated")
proj <- loadArchRProject(save_dir)
table(proj$Sample, proj$integrated_subclusters)

## add the pseudobulk peaks by cell type
proj <- addGroupCoverages(ArchRProj = proj, groupBy = "integrated_subclusters", 
                          minCells = 50, maxCells = 500, maxFragments = 25 * 10^6,
                          minReplicates = 2, maxReplicates = 12, force = T)
proj = saveArchRProject(ArchRProj = proj)

## call peaks w/ macs2, using reproducibility of peaks more than 
proj = proj %>% 
  addReproduciblePeakSet(groupBy = "integrated_subclusters", minCells = 50, force = T, 
  reproducibility = "(n+1)/4", peaksPerCell = 500, maxPeaks = 300000, plot = F)
proj = saveArchRProject(ArchRProj = proj)

## add peak matrix
proj <- addPeakMatrix(proj)
proj = saveArchRProject(ArchRProj = proj)


######################################################
## 3) get the reproducible peaks across cell types
source('/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate/code/raw_code/hal_scripts/narrowPeakFunctions.R')

LABEL='Human_Striatum_snATAC'; 
GENOME = 'hg38'; 

#######################################################
# get the reproducible peaks across Cluster2 cell types
peak_rds_fn = list.files(path = here(save_dir,'PeakCalls'), full.names = T, pattern = '.rds')
names(peak_rds_fn) = ss(basename(peak_rds_fn),'-reproduciblePeaks.gr.rds')
peakList = lapply(peak_rds_fn, readRDS)

# label the summit and peak name using hg38 coordinates
peakList = lapply(peakList, addSummitCenter)
peakList = lapply(peakList, nameNarrowPeakRanges, genome = GENOME)
peakList = lapply(peakList, sort)

lengths(peakList) # for processed dataset notebook

###############################################
# create directory and narrowPeak file names ##
PEAKDIR2= here('data/raw_data/peak_hg38/')
system(paste('mkdir -p',PEAKDIR2))
narrowPeak_hg38_fn = paste0(PEAKDIR2, paste0(LABEL, '.', names(peakList), 
                                             '.hg38.narrowPeak.gz'))
# write peaks to narrowPeak file
outRanges = mapply(write_GRangesToNarrowPeak,
                   gr = peakList, file = narrowPeak_hg38_fn, genome = GENOME)

###############################################
varFeats = c(300)

## clustering
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
      sampleCellsPre = 50000,
      varFeatures = varFeat*1000,
      saveIterations = FALSE,
      useMatrix = "PeakMatrix", 
      depthCol = "nFrags",
      name = iterLSIName,
      force = TRUE
    )

    proj <- saveArchRProject(proj)

    #batchcorrection
    HarmonyName = paste0("HarmonyX",varFeat)
    proj <- addHarmony(proj, reducedDims = iterLSIName,
                       name = HarmonyName, groupBy = 
                         c("Sample","Region", 'Publication', 'Donor'), 
                       lambda = rep(1, 4),
                       force = TRUE)
    #clustering
    ClustersName = paste0("ClustersX",varFeat)
    proj <- addClusters(proj, reducedDims = HarmonyName, name = ClustersName, 
                        resolution = 1, maxClusters = 35, force = TRUE)

    #Add UMAP
    UMAPName = paste0("UMAPH",varFeat)
    proj <- addUMAP(
      ArchRProj = proj, reducedDims = HarmonyName, name = UMAPName, 
      nNeighbors = 30, minDist = 0.5, metric = "cosine", force = TRUE)
    proj <- addImputeWeights(proj, reducedDims = HarmonyName)
    proj = saveArchRProject(proj)
}

