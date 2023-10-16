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


# Li et al.CATlas BICCN
mouse_fn2 = here(DATADIR,'Mouse_Li/Li2021_mouse_striatum_msns.hgGenes.h5Seurat')
obj_mm3 = LoadH5Seurat(mouse_fn2)

obj_mm3@meta.data = obj_mm3[[]] %>% 
  mutate(Region = case_when(
    grepl('CP', Sub.Region) ~ 'CPu',
    grepl('ACB', Sub.Region) ~ 'NAcc/OT'), 
    Species = "Mouse", 
    Project = 'Li_et_al',
    celltype.orig = SubClass,
    Replicate = orig.ident)

table(obj_mm3$celltype.orig)

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

# Corces et al. Caud/Put human CTL only dataset
corces_fn = here(DATADIR,'Human_Corces/Corces2020_human_striatum_msns.hgGenes.h5Seurat')
obj_corces = LoadH5Seurat(corces_fn)
obj_corces$Species = "Human"
obj_corces$Project = "Corces_et_al"
obj_corces$Region = "Caudate"
obj_corces$Replicate = obj_corces$Subject
obj_corces$celltype.orig = obj_corces$ClusterName

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
  ## integrate the human RNA then ATAC datasets together first
  'Human_Tran' = obj_hgNAc, 'Human_Phan' = obj_hgDStr, 'Human_Corces' = obj_corces,
  ## quality between the 2 monkeys so bad, need to integrate the 2 monkeys together
  'Macaque_HeKleyman1' = obj_rm1, 'Macaque_HeKleyman2' = obj_rm2, 
  ## then integrate mouse RNA and ATAC datasets together
  'Mouse_Stanley' = obj_mm, 'Mouse_Saunders' = obj_mm2, 'Mouse_Li' = obj_mm3, 
  ## then integrate with the Rat
  'Rat_Savell' = obj_rn)

## specify order of integration defined by
## https://satijalab.org/seurat/reference/integratedata
integrate.matrix = matrix(c(-1, -2, # the humans, RNA
                             1, -3, # the human RNA + ATAC
                            -4, -5, # the monkeys
                             2,  3, # the primates
                            -6, -7, # the mice, RNA
                             5, -8, # the mouse RNA + ATAC
                            -9,  6, # the rodents
                             4,  7  # the mammals
), ncol = 2, byrow = T)

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
  FindClusters(resolution = .8, algorithm = 2)

###############################################################
## 3) plot the seurat clusters to see which clusters are which

## relabel some things for organisation
proj_lvl = c('Tran_et_al', 'Phan_et_al', 'Corces_et_al', 'HeKleyman_et_al',
             'Stanley_et_al', 'Saunders_et_al', 'Li_et_al', 'Savell_et_al')
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
    width = 11, height = 5)
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


###############################################################
## 4) assign the consensus labels that are stable across species

## manually update clusters to be the main and sub clusters
combined.sct@meta.data = combined.sct[[]] %>% 
  mutate(
    Region = ifelse(Region == 'CAUD', 'Caudate', Region),
   integrated_clusters = case_when(seurat_clusters %in% c(0,2,4,5,11,15,18) ~ 'D1', 
                                   seurat_clusters %in% c(1,3,7,9,12,13,16,19) ~ 'D2',  
                                   seurat_clusters %in% c(6,8,10) ~ 'eSPN',  
                                   seurat_clusters %in% c(14) ~ 'IC', 
                                   seurat_clusters %in% c(17) ~ 'Other'), 
   
   integrated_subclusters = case_when(seurat_clusters %in% c(0,2,5,11,15,18) ~ 'D1-Matrix', 
                                      seurat_clusters %in% c(4) ~ 'D1-Striosome', 
                                      seurat_clusters %in% c(1,3,9,13,16,19) ~ 'D2-Matrix',  
                                      seurat_clusters %in% c(7,12) ~ 'D2-Striosome',  
                                      seurat_clusters %in% c(6) ~ 'D1/D2H', 
                                      seurat_clusters %in% c(8,10) ~ 'D1-NUDAP', 
                                      seurat_clusters %in% c(14) ~ 'IC', 
                                      seurat_clusters %in% c(17) ~ 'Other'))

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
cols = setNames(RColorBrewer::brewer.pal(8,'Paired'), 
                c('D1-Matrix', 'D1-Striosome', 'D1-NUDAP',  'D1/D2H', 
                  'D2-Matrix', 'D2-Striosome', 'IC', 'Other'))

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

## drop the low QC cluster
combined.sct = subset(combined.sct, integrated_clusters != 'Other')

#####################################################
## 5) save the integrated object w/ ATAC and RNA data
save_fn = here(DATADIR,'rdas/integrated_mouse_rat_macaque_human_striatum.ARC.h5Seurat')
SaveH5Seurat(combined.sct, save_fn, overwrite = T)

## save the integrated object with just the RNA
combined.sct2 = combined.sct[,!combined.sct$Project %in% c('Corces et al.', 'Li et al.')]
save_fn2 = here(DATADIR,'rdas/integrated_mouse_rat_macaque_human_striatum.RNA.h5Seurat')
SaveH5Seurat(combined.sct2, save_fn2, overwrite = T)

## save the integrated object with just the RNA
combined.sct3 = combined.sct[,combined.sct$Project %in% c('Corces et al.', 'Li et al.')]
save_fn3 = here(DATADIR,'rdas/integrated_mouse_rat_macaque_human_striatum.ATAC.h5Seurat')
SaveH5Seurat(combined.sct3, save_fn3, overwrite = T)



