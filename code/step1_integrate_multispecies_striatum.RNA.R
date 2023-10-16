library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(data.table)
library(future)
library(viridis)
library(harmony)
library(Nebulosa)
library(here)

# Enable parallelization
plan("multicore", workers = 6)
options(future.globals.maxSize = 100 * 1024^3, future.seed=TRUE)

DATADIR = 'data/tidy_data/'
PLOTDIR = 'figures/multispecies_integration/plots'
dir.create(PLOTDIR, recursive = T)

##############################
## 1) load in the datasets

# Stanley et al. mouse dataset w/ human genes
mouse_fn = here(DATADIR,'Mouse_Stanley/Stanley2020_mouse_striatum_msns.hgGenes.h5Seurat')
obj_mm = LoadH5Seurat(mouse_fn)
table(obj_mm$continuous.subtype, obj_mm$discrete.type)

obj_mm@meta.data = obj_mm[[]] %>% 
  mutate(Region = case_when(
    grepl('^Pat|^Mat|CPu|pat', continuous.subtype) ~ 'CPu',
    grepl('Shell|NAcc|ICj|^OT|vClust|unknown', continuous.subtype) ~ 'NAcc/OT'), 
    Species = "Mouse", Project = 'Stanley_et_al',
    celltype.orig = discrete.type, 
    Replicate = expt.name)

table(obj_mm$celltype.orig, obj_mm$Region)

# Saunders et al. Dropviz mouse dataset w/ human genes
mouse_fn2 = here(DATADIR,'Mouse_Saunders/Saunders2018_mouse_striatum_msns.hgGenes.h5Seurat')
obj_mm2 = LoadH5Seurat(mouse_fn2)

obj_mm2@meta.data = obj_mm2[[]] %>% 
  mutate(Region = 'CPu',
    Species = "Mouse", Project = 'Saunders_et_al',
    celltype.orig = str_replace_all(common_name, ' markers', '') %>% 
      str_replace_all(', Th+', '') %>% str_replace_all('\\+',''),
    Replicate = orig.ident)

table(obj_mm2$celltype.orig)

# rat dataset w/ human genes
rat_fn = here(DATADIR,'Rat_Savell/Savell2020_rat_striatum_msns.hgGenes.h5Seurat')
obj_rn = LoadH5Seurat(rat_fn)
table(obj_rn$CellType)

obj_rn@meta.data = obj_rn[[]] %>% 
  mutate(Region = 'NAcc', 
    Species = "Rat", Project = 'Savell_et_al',
    celltype.orig = CellType, 
    Replicate = Sample)

# Tran et al. NAcc human dataset
tran_fn = here(DATADIR,'Human_Tran/SCE_NAc-n8_tran-etal.msns.h5Seurat')
obj_hgNAc = LoadH5Seurat(tran_fn)
obj_hgNAc = RenameAssays(object = obj_hgNAc, originalexp = 'RNA')
obj_hgNAc$Region = "NAcc"
obj_hgNAc$Species = "Human"
obj_hgNAc$Project = "Tran_et_al"
obj_hgNAc$Replicate = obj_hgNAc$donor
obj_hgNAc$celltype.orig = as.character(obj_hgNAc$cellType)

# Phan et al. Caud/Put human CTL only dataset
phan_fn = here(DATADIR,'Human_Phan/OUD_Striatum_refined_SeuratObj_N10.msns.h5Seurat')
obj_hgDStr = LoadH5Seurat(phan_fn)
DefaultAssay(obj_hgDStr) = 'RNA'
obj_hgDStr$Species = "Human"
obj_hgDStr$Project = "Phan_et_al"
obj_hgDStr$Replicate = obj_hgDStr$ID
obj_hgDStr$celltype.orig = obj_hgDStr$celltype3

# He, Kleyman et al. Caud/Put/Nacc macaque dataset
macaque_fn = here(DATADIR,'Macaque_HeKleyman/GSE167920_Results_MSNs_processed_final.h5Seurat')
obj_rm = LoadH5Seurat(macaque_fn)
DefaultAssay(obj_rm) = 'RNA'
obj_rm@meta.data = obj_rm[[]] %>% 
  mutate(Region = case_when(
    region_name == 'caudate' ~ "Caudate", 
    region_name == 'putamen' ~ 'Putamen', 
    region_name == 'nacc' ~ 'NAcc'), 
    Species = "Macaque", Project = "HeKleyman_et_al",
    Replicate = monkey,
    celltype.orig = as.character(MSN_type))
obj_rm1 = subset(obj_rm, Replicate == 'Monkey_F')
obj_rm2 = subset(obj_rm, Replicate == 'Monkey_P')

#######################################################
## 2) integrate the datasets by shared common genes
obj_list = list(
  ## integrate the human datasets together first
  'Human_Tran' = obj_hgNAc, 'Human_Phan' = obj_hgDStr, 
  ## quality between the 2 monkeys so bad, need to integrate the 2 monkeys together
  'Macaque_HeKleyman1' = obj_rm1, 'Macaque_HeKleyman2' = obj_rm2, 
  ## then lastly integrate the rat and mouse datasets
  'Mouse_Stanley' = obj_mm, 'Mouse_Saunders' = obj_mm2, 'Rat_Savell' = obj_rn)

# subset to just the genes that overlap all datasets
genes = obj_list %>% lapply(rownames) %>% Reduce(f = 'intersect')
length(genes) #11740 genes

## create the multiple object list and perform SCTransform normaliation
obj_list = obj_list %>% lapply(function(x) x[genes,]) # subset to just shared genes
obj_list <- lapply(X = obj_list, FUN = SCTransform, method = "glmGamPoi",vst.flavor = "v2")
features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = 3000)
obj_list <- PrepSCTIntegration(object.list = obj_list, anchor.features = features)
obj_list <- lapply(X = obj_list, FUN = RunPCA, features = features, verbose = F)

## integrate normalized objects
anchors <- FindIntegrationAnchors(obj_list, normalization.method = "SCT",
                                  reduction = "rpca", k.anchor = 25,
                                  anchor.features = features, dims = 1:20)

## specify order of integration defined by
## https://satijalab.org/seurat/reference/integratedata
integrate.matrix = matrix(c(-1, -2, # the humans
                            -3, -4, # the monkeys
                             1,  2, # the primates
                            -5, -6, # the mice
                            -7,  4, # the rodents
                             3,  5  # the mammals
                            ), ncol = 2, byrow = T)

## do the integration of multiple datasets together
combined.sct = IntegrateData(anchorset = anchors, normalization.method = "SCT", 
                             dims = 1:20, sample.tree = integrate.matrix) %>% 
  RunPCA(verbose = F, assay = 'integrated')

# use Harmony to finish out the batch correction for different samples
combined.sct = combined.sct %>% 
  RunHarmony(c("Project", "Replicate"), max_iter = 20, lambda = c(1,1)) %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:30) %>% 
  # test this is the best resolution b/t lump & split
  FindClusters(resolution = 1, algorithm = 2)

###############################################################
## 3) plot the seurat clusters to see which clusters are which

## relabel some things for organisation
proj_lvl = c('Tran_et_al', 'Phan_et_al', 'HeKleyman_et_al', 'Stanley_et_al',
             'Saunders_et_al', 'Savell_et_al')
combined.sct$Project = factor(combined.sct$Project, proj_lvl)
levels(combined.sct@meta.data$Project) = levels(combined.sct$Project) %>% 
  str_replace_all('_', ' ') %>% str_replace_all('HeK', 'He, K') %>% 
  paste0('.')

rep_lvl = combined.sct[[]] %>% arrange(Project, Replicate) %>% 
  pull(Replicate) %>% unique()
combined.sct$Replicate = factor(combined.sct$Replicate, rep_lvl)

cell_lvl = combined.sct[[]] %>% arrange(Project, celltype.orig) %>% 
  pull(celltype.orig) %>% unique()
combined.sct$celltype.orig = factor(combined.sct$celltype.orig, cell_lvl)

## make the plots
pdf(here::here(PLOTDIR, 'Replicates_plots_bySpecies_clusters.umap.pdf'), 
    width = 8, height = 5)
DimPlot(combined.sct, reduction = "umap", split.by = "Species", 
        cols = ArchR::paletteDiscrete(combined.sct$seurat_clusters),
        label = T,
        group.by = 'seurat_clusters') + theme(legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 5, override.aes = list(size = 3))) &
  FontSize(main = 7, x.title = 7, y.title = 7)
dev.off()

pdf(here::here(PLOTDIR, 'Replicates_plots_byProject_origLabel.umap.pdf'), 
    width = 8, height = 5)
DimPlot(combined.sct, reduction = "umap", split.by = "Project", 
        cols = ArchR::paletteDiscrete(combined.sct$celltype.orig),
        label = T,
        group.by = 'celltype.orig') + theme(legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 6, override.aes = list(size = 3))) &
  FontSize(main = 7, x.title = 7, y.title = 7)
dev.off()

## manually update clusters to be the main and sub clusters
combined.sct@meta.data = combined.sct[[]] %>% 
  mutate(
   integrated_clusters = case_when(seurat_clusters %in% c(1,2,5,7,10,12,15,18) ~ 'D2',  
                                   seurat_clusters %in% c(0,3,4,6,9,11,16,17) ~ 'D1', 
                                   seurat_clusters %in% c(8,13) ~ 'eSPN',  
                                   seurat_clusters %in% c(14) ~ 'IC'), 
   
   integrated_subclusters = case_when(seurat_clusters %in% c(1,2,5,10,12,18) ~ 'D2-Matrix',  
                                      seurat_clusters %in% c(7,15) ~ 'D2-Striosome',  
                                      seurat_clusters %in% c(0,3,4,9,11,16,17) ~ 'D1-Matrix', 
                                      seurat_clusters %in% c(6) ~ 'D1-Striosome', 
                                      seurat_clusters %in% c(8) ~ 'D1/D2H', 
                                      seurat_clusters %in% c(13) ~ 'D1-NUDAP', 
                                      seurat_clusters %in% c(14) ~ 'IC'))

table(is.na(combined.sct$integrated_clusters), combined.sct$Project)
table(combined.sct$seurat_clusters, combined.sct$integrated_clusters)

## make the plots
pdf(here::here(PLOTDIR, 'Replicates_plots_bySpecies_mainCluster.umap.pdf'), 
    width = 8, height = 5)
DimPlot(combined.sct, reduction = "umap", split.by = "Species", 
        cols = ArchR::paletteDiscrete(combined.sct$integrated_clusters),
        label = T,
        group.by = 'integrated_clusters') + theme(legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 5, override.aes = list(size = 3))) &
  FontSize(main = 7, x.title = 7, y.title = 7)
dev.off()

## make the plots
cols = setNames(RColorBrewer::brewer.pal(7,'Paired'), 
                c('D1-Matrix', 'D1-Striosome', 'D1-NUDAP',  'D1/D2H', 
                  'D2-Matrix', 'D2-Striosome', 'IC'))

combined.sct@meta.data$integrated_subclusters = 
  factor(combined.sct$integrated_subclusters, names(cols))

pdf(here::here(PLOTDIR, 'Replicates_plots_bySpecies_subCluster.umap.pdf'), 
    width = 8, height = 5)
DimPlot(combined.sct, reduction = "umap", split.by = "Species", 
        cols = cols, label = T, group.by = 'integrated_subclusters') +
  theme(legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 5, override.aes = list(size = 3))) &
  FontSize(main = 7, x.title = 7, y.title = 7)
dev.off()

##################################
## 4) save the integrated object
save_fn = here(DATADIR,'rdas/integrated_mouse_rat_macaque_human_striatum.h5Seurat')
SaveH5Seurat(combined.sct, save_fn, overwrite = T)

