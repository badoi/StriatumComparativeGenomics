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

# 1, Stanley et al. mouse dataset w/ human genes
mouse_fn = here(DATADIR,'Mouse_Stanley/Stanley2020_mouse_striatum_msns.hgGenes.h5Seurat')
obj_mm1 = LoadH5Seurat(mouse_fn)
table(obj_mm1$continuous.subtype, obj_mm1$discrete.type)

obj_mm1@meta.data = obj_mm1[[]] %>% 
  mutate(Region = case_when(
    grepl('^Pat|^Mat|CPu|pat', continuous.subtype) ~ 'CPu',
    grepl('Shell|NAcc|ICj|^OT|vClust|unknown', continuous.subtype) ~ 'NAcc/OT'), 
    Species = "Mouse", Project = 'Stanley_et_al',
    Publication.Celltype = discrete.type, 
    Replicate = expt.name, 
    Case = 'Unaffected', Sex = 'Male')

table(obj_mm1$Publication.Celltype, obj_mm1$Region)

# 2, Saunders et al. Dropviz mouse dataset w/ human genes
mouse_fn2 = here(DATADIR,'Mouse_Saunders/Saunders2018_mouse_striatum_msns.hgGenes.h5Seurat')
obj_mm2 = LoadH5Seurat(mouse_fn2)

obj_mm2@meta.data = obj_mm2[[]] %>% 
  mutate(Region = 'CPu',
         Species = "Mouse", 
         Project = 'Saunders_et_al',
         Publication.Celltype = str_replace_all(common_name, ' markers', '') %>% 
           str_replace_all(', Th+', '') %>% str_replace_all('\\+',''),
         Replicate = orig.ident, 
         Case = 'Unaffected', 
         Sex = 'Male')

table(obj_mm2$Publication.Celltype)

# 3, rat dataset w/ human genes, acute & repeated cocaine exposure
obj_rn3 = here(DATADIR,'Rat_Phillips', 
               'Phillips2023_snRNA_filtered_SeuratObj_N8.hgGenes.rds') %>% readRDS()

obj_rn3@meta.data = obj_rn3[[]] %>% 
  mutate(Region = 'NAcc', 
         Species = "Rat", 
         Project = ifelse(Dataset == 'Repeated', 'Phillips_et_al', 'Savell_et_al'),
         Publication.Celltype = Combo_CellType, 
         Replicate = Sample, 
         Case = ifelse(Stim == 'Coc' & Dataset == 'Repeated', 'Cocaine-Repeated', 
                       ifelse(Stim == 'Coc' & Dataset == 'Acute', 'Cocaine-Acute','Saline')), 
         Sex = ifelse(Sex == 'Fem', 'Female', 'Male'))

table(obj_rn3$Publication.Celltype)
table(obj_rn3$Case)
table(obj_rn3$Sex)

# 4,5, He, Kleyman et al. Caud/Put/Nacc macaque dataset
macaque_fn = here(DATADIR,'Macaque_HeKleyman/GSE167920_Results_MSNs_processed_final.h5Seurat')
obj_rm = LoadH5Seurat(macaque_fn, assay = 'RNA')
obj_rm@meta.data = obj_rm[[]] %>% 
  mutate(Region = case_when(
    region_name == 'caudate' ~ "Caudate", 
    region_name == 'putamen' ~ 'Putamen', 
    region_name == 'nacc' ~ 'NAcc'), 
    Species = "Macaque", 
    Project = "HeKleyman_et_al",
    Replicate = monkey,
    Sex = ifelse(Replicate == 'Monkey_F','Female', 'Male'),
    Publication.Celltype = as.character(MSN_type), 
    Case = 'Unaffected')

# these two have such different qualities, need to be treated separately
obj_rm4 = subset(obj_rm, Replicate == 'Monkey_F')  
obj_rm5 = subset(obj_rm, Replicate == 'Monkey_P')

# 6, Phan et al. Caud/Put human dataset
phan_fn = here(DATADIR,'Human_Phan/BU_OUD_Striatum_refined_msn_SeuratObj_N22.h5Seurat')
obj_hg6 = LoadH5Seurat(phan_fn, assay = 'RNA') %>% 
  FindVariableFeatures(verbose = F) %>% ScaleData(verbose = F) %>% RunPCA(verbose = F) 
obj_hg6$Species = "Human"
obj_hg6$Project = "Phan_et_al"
obj_hg6$Replicate = obj_hg6$ID
obj_hg6$Publication.Celltype = obj_hg6$celltype3 %>% as.character()
obj_hg6$Case = ifelse(obj_hg6$DSM.IV.OUD == 'CON', 'Unaffected', 'Opioid use disorder')
obj_hg6$Sex = ifelse(obj_hg6$Sex == 'M', 'Male', 'Female')
table(obj_hg6$Region)

# 7, Tran et al. NAcc human dataset
tran_fn = here(DATADIR,'Human_Tran/SCE_NAc-n8_tran-etal.msns.h5Seurat')
obj_hg7 = LoadH5Seurat(tran_fn) %>% 
  FindVariableFeatures(verbose = F) %>% ScaleData(verbose = F) %>% RunPCA(verbose = F) 
obj_hg7 = RenameAssays(object = obj_hg7, originalexp = 'RNA')
obj_hg7$Region = "NAcc"
obj_hg7$Species = "Human"
obj_hg7$Project = "Tran_et_al"
obj_hg7$Replicate = obj_hg7$donor
obj_hg7$Publication.Celltype = obj_hg7$cellType %>% as.character()
obj_hg7$Case = 'Unaffected'
obj_hg7$Sex = ifelse(obj_hg7$sex == 'M', 'Male', 'Female')

#######################################################
## 2) integrate the datasets by shared common genes
obj_list = list(
  ## then integrate mouse RNA datasets together
  'Mouse_Stanley' = obj_mm1, 'Mouse_Saunders' = obj_mm2,  
  ## then integrate with the Rat
  'Rat_Phillips' = obj_rn3,
  ## quality between the 2 monkeys so bad, need to split the monkeys apart
  'Macaque_HeKleyman1' = obj_rm4, 'Macaque_HeKleyman2' = obj_rm5,
  ## integrate the human RNA datasets together first
  'Human_Phan' = obj_hg6, 'Human_Tran' = obj_hg7
)

## specify order of integration defined by
## https://satijalab.org/seurat/reference/integratedata
integrate.matrix = matrix(c(-1, -2, # 1: the mice
                            -3, 1,  # 2: the rodents
                            -4, -5, # 3: the monkeys
                            -6, -7, # 4: the humans
                            3,  4,  # 5: the primates, outputs of lines 3 & 4
                            2,  5   # 6: the mammals, outputs of lines 2 & 5
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
                                  reduction = "rpca", k.anchor = 40,
                                  anchor.features = features, dims = 1:20)

## do the integration of multiple datasets together
combined.sct = IntegrateData(anchorset = anchors, normalization.method = "SCT", 
                             dims = 1:20, sample.tree = integrate.matrix)

# use Harmony to finish out the batch correction for different samples
combined.sct = combined.sct %>% 
  RunPCA(verbose = F, assay = 'integrated') %>% 
  RunHarmony(c("Project", "Replicate", 'Species', 'Sex', 'Case', 'Region'), 
             max_iter = 20, lambda = rep(1,6)) %>% 
  RunUMAP(reduction = "harmony", dims = 1:30, return.model = T) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:30) %>% 
  # test this is the best resolution b/t lump & split
  FindClusters(resolution = 1, algorithm = 2)

###############################################################
## 3) plot the seurat clusters to see which clusters are which

## relabel some things for organisation
proj_lvl = c('Tran_et_al', 'Phan_et_al', 'HeKleyman_et_al', 'Stanley_et_al',
             'Saunders_et_al', 'Phillips_et_al', 'Savell_et_al')
combined.sct$Project = factor(combined.sct$Project, proj_lvl)

levels(combined.sct@meta.data$Project) = levels(combined.sct$Project) %>%
  str_replace_all('_', ' ') %>% str_replace_all('HeK', 'He, K') %>%
  paste0('.')

rep_lvl = combined.sct[[]] %>% arrange(Project, Replicate) %>% 
  pull(Replicate) %>% unique()
combined.sct$Replicate = factor(combined.sct$Replicate, rep_lvl)

cell_lvl = combined.sct[[]] %>% arrange(Project, Publication.Celltype) %>% 
  pull(Publication.Celltype) %>% unique()
combined.sct$Publication.Celltype = factor(combined.sct$Publication.Celltype, cell_lvl)

## make the plots
pdf(here::here(PLOTDIR, 'Replicates_plots_bySpecies_clusters.umap.pdf'), 
    width = 11, height = 5)
DimPlot(combined.sct, reduction = "umap", split.by = "Species", 
        cols = ArchR::paletteDiscrete(combined.sct$seurat_clusters),
        label = T,
        group.by = 'seurat_clusters') + theme(legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 2, override.aes = list(size = 3))) &
  FontSize(main = 7, x.title = 7, y.title = 7)
dev.off()

pdf(here::here(PLOTDIR, 'Replicates_plots_byProject_origLabel.umap.pdf'), 
    width = 11, height = 5)
DimPlot(combined.sct, reduction = "umap", split.by = "Project", 
        cols = ArchR::paletteDiscrete(combined.sct$Publication.Celltype),
        label = T,
        group.by = 'Publication.Celltype') + theme(legend.position = 'bottom') + 
  guides(colour = guide_legend(nrow = 6, override.aes = list(size = 3))) &
  FontSize(main = 7, x.title = 7, y.title = 7)
dev.off()

## manually update clusters to be the main and sub clusters
combined.sct@meta.data = combined.sct[[]] %>% 
  mutate(
   integrated_clusters = case_when(
     seurat_clusters %in% c(1,4,5,8,10,12,14,16,17,19,20) ~ 'D2',  
     seurat_clusters %in% c(0,2,3,6,9,15,18,21) ~ 'D1', 
     seurat_clusters %in% c(7,11) ~ 'eSPN', 
     seurat_clusters %in% c(13) ~ 'IC'), 
   
   integrated_subclusters = case_when(
     seurat_clusters %in% c(1,4,5,8,16,19,20) ~ 'D2-Matrix',  
     seurat_clusters %in% c(10,12,14,17) ~ 'D2-Striosome',  
     seurat_clusters %in% c(0,2,3,9,15,21) ~ 'D1-Matrix', 
     seurat_clusters %in% c(6,18) ~ 'D1-Striosome', 
     seurat_clusters %in% c(7) ~ 'D1/D2H', 
     seurat_clusters %in% c(11) ~ 'D1-NUDAP', 
     seurat_clusters %in% c(13) ~ 'IC')
   )

table(is.na(combined.sct$integrated_clusters), combined.sct$Project)
table(is.na(combined.sct$integrated_clusters), combined.sct$seurat_clusters)
table(is.na(combined.sct$integrated_subclusters), combined.sct$seurat_clusters)

## make the plots
pdf(here::here(PLOTDIR, 'Replicates_plots_bySpecies_mainCluster.umap.pdf'), 
    width = 8, height = 5)
DimPlot(combined.sct, reduction = "umap", split.by = "Species", 
        cols = ArchR::paletteDiscrete(combined.sct$integrated_clusters),
        label = T, group.by = 'integrated_clusters') + 
  theme(legend.position = 'bottom') + 
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

## once this is correct, then drop the NA clusters
combined.sct = combined.sct[, !is.na(combined.sct$integrated_clusters)]

##################################
## 4) save the integrated object
save_fn = here(DATADIR,'rdas/integrated_mouse_rat_macaque_human_striatum.RNA.rds')
SaveSeuratRds(combined.sct, save_fn)


