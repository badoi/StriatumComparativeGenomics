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
TEMPDIR = '/scratch/bnphan'
PLOTDIR = 'figures/multispecies_integration/plots'
dir.create(PLOTDIR, recursive = T)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

##############################
## 1) load in the datasets

# Stanley et al. mouse dataset w/ human genes
mouse_fn = here(DATADIR,'Mouse_Stanley/Stanley2020_mouse_striatum_msns.hgGenes.h5Seurat')
obj_mm1 = LoadH5Seurat(mouse_fn)
table(obj_mm1$continuous.subtype, obj_mm1$discrete.type)

obj_mm1@meta.data = obj_mm1[[]] %>% 
  mutate(Region = case_when(
    grepl('^Pat|^Mat|CPu|pat', continuous.subtype) ~ 'CPu',
    grepl('Shell|NAcc|ICj|^OT|vClust|unknown', continuous.subtype) ~ 'NAcc'), 
    Species = "Mouse", 
    Assay = "RNA",
    Project = 'Stanley et al.',
    Publication.Celltype = discrete.type, 
    Replicate = expt.name, 
    Case = 'Unaffected', Sex = 'Male')

table(obj_mm1$Publication.Celltype, obj_mm1$Region)

# Saunders et al. Dropviz mouse dataset w/ human genes
mouse_fn2 = here(DATADIR,'Mouse_Saunders/Saunders2018_mouse_striatum_msns.hgGenes.h5Seurat')
obj_mm2 = LoadH5Seurat(mouse_fn2)

obj_mm2@meta.data = obj_mm2[[]] %>% 
  mutate(Region = 'CPu',
         Species = "Mouse", 
         Project = 'Saunders et al.',
         Assay = "RNA",
         Publication.Celltype = str_replace_all(common_name, ' markers', '') %>% 
           str_replace_all(', Th+', '') %>% str_replace_all('\\+',''),
         Replicate = orig.ident, 
         Case = 'Unaffected', 
         Sex = 'Male')

table(obj_mm2$Publication.Celltype)

# Chen et al. Nat Neuro NAcc mouse dataset w/ human genes
obj_mm3 = here(DATADIR,'Mouse_Chen/Chen2021_mouse_NAc_snRNA_msn.hgGenes.rds') %>% 
  readRDS()

obj_mm3@meta.data = obj_mm3[[]] %>% 
  mutate(Region = 'NAcc',
         Species = "Mouse", 
         Project = 'Chen et al.',
         Assay = "RNA",
         Publication.Celltype = CellType,
         Replicate = Sample, 
         Case = 'Unaffected', 
         Sex = 'Male')

table(obj_mm3$Publication.Celltype)

# BICCN Allen brain subset 2024
obj_mm4 = here(DATADIR,'Mouse_Zeng/Zeng2023_mouse_striatum_msns.hgGenes.rds') %>% 
  readRDS() 

obj_mm4@meta.data = obj_mm4[[]] %>% 
  mutate(Region = case_when(region_of_interest_acronym == 'STRd' ~ 'CPu', 
                            region_of_interest_acronym == 'STRv' ~ 'NAcc'),
         Species = "Mouse", 
         Project = 'Zeng et al.',
         Assay = "RNA",
         Publication.Celltype = subclass %>% str_remove('([0-9]+ )'),
         Replicate = donor_label, 
         Case = 'Unaffected', 
         Sex = ifelse(donor_sex == 'F','Female' , 'Male'))

table(obj_mm4$Publication.Celltype)

#####################################################################
## rat dataset w/ human genes, acute & repeated cocaine exposure
obj_rn1 = here(DATADIR,'Rat_Phillips', 
               'Phillips2023_snRNA_filtered_SeuratObj_N8.hgGenes.rds') %>% readRDS()

obj_rn1 = obj_rn1[,!is.na(obj_rn1$Combo_CellType)] 

obj_rn1@meta.data = obj_rn1[[]] %>% 
  mutate(Region = 'NAcc', 
         Species = "Rat", 
         Assay = "RNA",
         Project = ifelse(Dataset == 'Repeated', 'Phillips et al.', 'Savell et al.'),
         Publication.Celltype = Combo_CellType, 
         Replicate = Sample, 
         Case = ifelse(Stim == 'Coc' & Dataset == 'Repeated', 'Cocaine-Repeated', 
                       ifelse(Stim == 'Coc' & Dataset == 'Acute', 'Cocaine-Acute','Saline')), 
         Sex = ifelse(Sex == 'Fem', 'Female', 'Male'))

table(obj_rn1$Publication.Celltype)
table(obj_rn1$Case)
table(obj_rn1$Sex)

#####################################################################
# He, Kleyman et al. Caud/Put/Nacc macaque dataset
macaque_fn = here(DATADIR,'Macaque_HeKleyman/GSE167920_Results_MSNs_processed_final.h5Seurat')
obj_rm = LoadH5Seurat(macaque_fn, assay = 'RNA')
obj_rm@meta.data = obj_rm[[]] %>% 
  mutate(Region = case_when( region_name == 'caudate' ~ "Caudate",  
                             region_name == 'putamen' ~ 'Putamen', 
                             region_name == 'nacc' ~ 'NAcc'), 
    Species = "Macaque", 
    Project = "HeKleyman et al.",
    Assay = "RNA",
    Replicate = monkey,
    Sex = ifelse(Replicate == 'Monkey_F','Female', 'Male'),
    Publication.Celltype = as.character(MSN_type), 
    Case = 'Unaffected')

# these two have such different qualities, need to be treated separately
obj_rm1 = subset(obj_rm, Replicate == 'Monkey_F')  
obj_rm2 = subset(obj_rm, Replicate == 'Monkey_P')

# Chiou et al. macaque striatum dataset
macaque_fn2 = here(DATADIR,'Macaque_Chiou/Chiou2023_macaque_sciRNA_basal_ganglia.msn.rds')
obj_rm3 = readRDS(macaque_fn2)

obj_rm3@meta.data = obj_rm3[[]] %>% 
  mutate(Region = case_when( region == 'CN' ~ "Caudate",  
                             region == 'NAc' ~ 'NAcc'), 
         Species = "Macaque", 
         Project = "Chiou et al.",
         Assay = "RNA",
         Replicate = as.character(donor_id),
         Sex = ifelse(sex == 'male', 'Male', 'Female'),
         Publication.Celltype = as.character(cell_subcluster), 
         Case = 'Unaffected')

table(obj_rm3$Publication.Celltype)
table(obj_rm3$Replicate)
table(obj_rm3$Sex)

##############################################
# Phan et al. Caud/Put human dataset
phan_fn = here(DATADIR,'Human_Phan/BU_OUD_Striatum_refined_msn_SeuratObj_N22.h5Seurat')
obj_hg1 = LoadH5Seurat(phan_fn, assay = 'RNA') %>% 
  FindVariableFeatures(verbose = F) %>% ScaleData(verbose = F) %>% RunPCA(verbose = F) 
obj_hg1$Species = "Human"
obj_hg1$Project = "Phan et al."
obj_hg1$Assay = "RNA"
obj_hg1$Replicate = obj_hg1$ID
obj_hg1$Publication.Celltype = obj_hg1$celltype3 %>% as.character()
obj_hg1$Case = ifelse(obj_hg1$DSM.IV.OUD == 'CON', 'Unaffected', 'Opioid use disorder')
obj_hg1$Sex = ifelse(obj_hg1$Sex == 'M', 'Male', 'Female')
table(obj_hg1$Region)

# Tran et al. NAcc human dataset
tran_fn = here(DATADIR,'Human_Tran/SCE_NAc-n8_tran-etal.msns.h5Seurat')
obj_hg2 = LoadH5Seurat(tran_fn) %>% 
  FindVariableFeatures(verbose = F) %>% ScaleData(verbose = F) %>% RunPCA(verbose = F) 
obj_hg2 = RenameAssays(object = obj_hg2, originalexp = 'RNA')
obj_hg2$Region = "NAcc"
obj_hg2$Species = "Human"
obj_hg2$Assay = "RNA"
obj_hg2$Project = "Tran et al."
obj_hg2$Replicate = obj_hg2$donor
obj_hg2$Publication.Celltype = obj_hg2$cellType %>% as.character()
obj_hg2$Case = 'Unaffected'
obj_hg2$Sex = ifelse(obj_hg2$sex == 'M', 'Male', 'Female')

# Silettie et al. BICCN human dataset
silettie_fn = here(DATADIR,'Human_Siletti/WHB-10Xv3-Neurons-Basal.nuclei.msn.seurat.rds')
obj_hg3 = readRDS(silettie_fn)

obj_hg3@meta.data = obj_hg3[[]] %>% 
  mutate(Region = case_when( region_of_interest_label == 'Human Pu' ~ 'Putamen', 
                             region_of_interest_label == 'Human CaB' ~ "Caudate",  
                             region_of_interest_label == 'Human NAC' ~ 'NAcc'),
         Species = "Human",
         Assay = "RNA",
         Project = "Siletti et al.",
         Replicate = donor_label,
         Publication.Celltype = name_supercluster,
         Case = 'Unaffected',
         Sex = ifelse(donor_sex == 'M', 'Male', 'Female'))

table(obj_hg3$Publication.Celltype)
table(obj_hg3$Replicate)
table(obj_hg3$Sex)

# Gayden et al. NAcc human dataset
gayden_fn = here(DATADIR,'Human_Gayden/Gayden_NAc_refined_msn_SeuratObj_N4.rds')
obj_hg6 = readRDS(gayden_fn)

obj_hg6$Region = "NAcc"
obj_hg6$Species = "Human"
obj_hg6$Assay = "RNA"
obj_hg6$Project = "Gayden et al."
obj_hg6$Replicate = obj_hg6$Case
obj_hg6$Publication.Celltype = obj_hg6$celltype3
obj_hg6$Case = 'Unaffected'
obj_hg6$Sex = ifelse(obj_hg6$Sex == 'M', 'Male', 'Female')
table(obj_hg6$Publication.Celltype)

#######################
## the ATAC datasets ##

# Zu et al.CATlas BICCN mouse
obj_mm5 = here(DATADIR,'Mouse_Zu/Zu2023_mouse_striatum_msns.hgGenes.rds') %>% 
  readRDS()

obj_mm5@meta.data = obj_mm5[[]] %>% 
  mutate(Region = case_when(
    grepl('CP', SubRegion) ~ 'CPu',
    grepl('ACB', SubRegion) ~ 'NAcc'), 
    Species = "Mouse", 
    Project = 'Zu et al.',
    Assay = 'ATAC',
    Publication.Celltype = 'InhibitoryNeurons',
    Replicate = Sample, 
    Sex = 'Male',
    Case = 'Unaffected')

table(obj_mm5$Publication.Celltype)

# Corces et al. human caudate
# Donor information from https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-020-00721-x/MediaObjects/41588_2020_721_MOESM3_ESM.xlsx
corces_fn = here(DATADIR,'Human_Corces/Corces2020_human_striatum_msns.hgGenes.h5Seurat')
obj_hg4 = LoadH5Seurat(corces_fn)
obj_hg4$Species = "Human"
obj_hg4$Project = "Corces et al."
obj_hg4$Assay = "ATAC"
obj_hg4$Region = "Caudate"
obj_hg4$Replicate = obj_hg4$Subject
obj_hg4$Publication.Celltype = obj_hg4$ClusterName
obj_hg4$Case = 'Unaffected'
obj_hg4$Sex = ifelse(obj_hg4$Replicate %in% c('09_1589', '14_1018'), 'Female', 'Male')

table(obj_hg4$Sex, obj_hg4$Replicate)

# Human CATLas BICCN dataset, Li et al. 2023
obj_hg5 = here(DATADIR,'Human_Li/Li2023_human_striatum_msns.rds') %>% readRDS()

obj_hg5@meta.data = obj_hg5[[]] %>% 
  mutate(Species = "Human",
         Project = "Li et al.",
         Assay = "ATAC",
         Region = case_when( Brain.region == 'Pu' ~ 'Putamen', 
                             Brain.region == 'CaB' ~ "Caudate",  
                             Brain.region == 'NAC' ~ 'NAcc'),
         Replicate = Brain.dissetion.ID %>% str_remove('\\.CX.*'),
         Publication.Celltype = case_when(grepl('D1|D2|FOX|MSN', celltype) ~ celltype, 
                                          TRUE ~ 'Other'),
         Publication.Celltype = ss(Publication.Celltype, '_', 1),
         Case = 'Unaffected',
         Sex = ifelse(grepl(c('H19.30.001|H19.30.002|H19.30.004'), Brain.dissetion.ID), 
                      'Male', 'Female')
         )
table(obj_hg5$Publication.Celltype, obj_hg5$Replicate)

#######################################################
## 2) integrate the datasets by shared common genes
obj_list = list(
  ## then integrate mouse RNA datasets together
  'Mouse_Stanley' = obj_mm1, # 1
  'Mouse_Saunders' = obj_mm2, # 2
  'Mouse_Chen' = obj_mm3, # 3
  ## the new baseline mouse dataset from Zeng et al. Allen Institute
  'Mouse_Zeng' = obj_mm4, # 4
  ## then integrate with the Rat
  'Rat_Phillips' = obj_rn1, # 5
  ## quality between the 2 monkeys so bad, need to split the monkeys apart
  'Macaque_HeKleyman1' = obj_rm1, # 6
  'Macaque_HeKleyman2' = obj_rm2, # 7
  ## then integrate the monkeys together
  'Macaque_Chiou' = obj_rm3, # 8
  ## integrate the human RNA datasets together first
  'Human_Phan' = obj_hg1, # 9
  'Human_Tran' = obj_hg2, # 10
  'Human_Gayden' = obj_hg6, # 11
  ## the new human dataset from Siletti et al. BICCN
  'Human_Siletti' = obj_hg3, # 12
  ## then integrate with the human ATAC datasets
  'Mouse_Zu' = obj_mm5, # 13
  'Human_Corces' = obj_hg4, # 14
  'Human_Li' = obj_hg5 # 15
)

# the names of obj_list in a vector written out by hand
obj_list_names = c('Mouse_Stanley', 'Mouse_Saunders', 'Mouse_Chen', 'Mouse_Zeng', 
                   'Rat_Phillips', 'Macaque_HeKleyman1', 'Macaque_HeKleyman2', 
                   'Macaque_Chiou', 'Human_Phan', 'Human_Tran', 'Human_Siletti', 
                   'Mouse_Zu', 'Human_Corces', 'Human_Li')
all.equal(names(obj_list), obj_list_names)

# subset to just the genes that overlap all datasets
genes = obj_list %>% lapply(rownames) %>% Reduce(f = 'intersect')
length(genes) #10200 genes shared across all datasets

## save the log normalized only files
files_list = names(obj_list) %>% paste0('.LogNorm.rds') %>% here(DATADIR, 'rdas/tmp', .)
files_list2 = names(obj_list) %>% paste0('.LogNorm.rds') %>% here(TEMPDIR, .)

obj_list2 = lapply(obj_list, function(obj){
  DefaultAssay(obj) = 'RNA'
  obj = obj %>% FindVariableFeatures(verbose = F) %>% 
    ScaleData(verbose = F) %>% RunPCA(verbose = F)
  obj[['SCT']] <- NULL
  return(obj)
})

parallel::mclapply(X = seq_along(obj_list2), function(i){
  saveRDS(obj_list2[[i]], files_list2[i])}, mc.cores = 14)

thecall = paste('rsync -Paq ', files_list2, ' ', files_list)
parallel::mclapply(thecall, system, mc.cores = 14)

#######################################################
## save this step because it takes a long time to run 
files_list = names(obj_list) %>% paste0('.sctransformed.rds') %>% here(DATADIR, 'rdas/tmp', .)
files_list2 = names(obj_list) %>% paste0('.sctransformed.rds') %>% here(TEMPDIR, .)
dir.create(dirname(files_list) %>% unique(), recursive = T, showWarnings = F)

## rewrite the following to run in series in a for loop
for( i in seq_along(obj_list)){
  if(!file.exists(files_list[i])){
    obj_list[[i]] = JoinLayers(obj_list[[i]])
    obj_list[[i]] <- SCTransform(obj_list[[i]], conserve.memory = T, verbose = T)
    saveRDS(obj_list[[i]], files_list2[i])
    system(paste('rsync -Paq ', files_list2[i], ' ', files_list[i]))
  } else {
    print(paste('File already exists:', files_list[i]))
  }
}
