library(tidyverse)
library(data.table)
library(SeuratDisk)
library(SuperCell)
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
cols1 = setNames(RColorBrewer::brewer.pal(7,'Paired'), subclusters)

# load the metadata of all cells
metadata = here('data/tidy_data/rdas/integrated_multispecies_striatum.meta.rds') %>% readRDS()
meta.rep = metadata %>% distinct(Project, Replicate, Assay, Species, Sex, Case, Region)
gamma <- 50 # Graining level

###################################################################
## 1) read in datasets, add the metadata, and create Supercells
data = c('Mouse_Chen/Chen2021_mouse_NAc_snRNA_msn.seurat.rds', 
         'Mouse_Saunders/Saunders2018_mouse_striatum_all.rds', 
         'Mouse_Stanley/Stanley2020_mouse_striatum_msns.rds', 
         'Mouse_Zeng/WMB-10Xv3-STR.msn.seurat.rds',
         'Rat_Phillips/Phillips2023_snRNA_msns_N8.rds', 
         
         'Human_Siletti/WHB-10Xv3-Neurons-Basal.nuclei.msn.seurat.rds', 
         'Human_Gayden/Gayden_NAc_refined_msn_SeuratObj_N4.rds', 
         'Human_Phan/BU_OUD_Striatum_refined_msn_SeuratObj_N22.h5Seurat',
         'Human_Tran/SCE_NAc-n8_tran-etal.msns.h5Seurat', 
         'Macaque_Chiou/Chiou2023_macaque_sciRNA_basal_ganglia.msn.rds', 
         'Macaque_HeKleyman/GSE167920_Results_MSNs_processed_final.rds')

outs = basename(data) %>% str_remove('.rds') %>% str_remove('.h5Seurat') %>% 
  paste0('.supercells.rds') %>% file.path(dirname(data), .)
  
for(i in seq_along(data)){
  # if(!file.exists(here(DATADIR, outs[i]))){
    x = data[i]
    # load the object 
    if(grepl('h5Seurat', x)) { obj = LoadH5Seurat(here(DATADIR, x))
    } else { obj = readRDS(here(DATADIR, x)) }
    
    # set the default assay to get counts
    if('RNA' %in% Assays(obj)) { DefaultAssay(obj) = 'RNA'
    } else if ('originalexp' %in% Assays(obj)) { DefaultAssay(obj) = 'originalexp'}
    
    # check if the counts are split
    if(!'counts' %in% Layers(obj)) { obj = JoinLayers(obj) }
    
    obj = obj[,Cells(obj) %in% metadata$CellID ]
    metadata2 = metadata[match(Cells(obj), metadata$CellID),]
    rownames(metadata2) = metadata2$CellID
    obj = AddMetaData(obj, metadata2)
    
    table(obj$Replicate, obj$integrated_subclusters)
    
    counts = LayerData(obj, "counts")
    rownames(counts) = str_replace(rownames(counts), 'Drd1a', 'Drd1')
    
    obj = CreateSeuratObject(counts, meta.data =obj[[]]) %>%
      NormalizeData() %>% FindVariableFeatures()

    # Compute metacells using SuperCell package
    MC <- SCimplify( X = LayerData(obj, "data"), genes.use = VariableFeatures(obj), 
                     cell.split.condition = obj$Replicate, gamma = gamma, n.pc = 5)
    
    MC.counts <- supercell_GE(ge = LayerData(obj, "counts"), mode = "sum", 
                              groups = MC$membership)
    MC.ge <- Seurat::LogNormalize(MC.counts, verbose = F)
    
    # transfer the metadata from the original object to the metacells
    MC$Replicate <- 
      supercell_assign(cluster = obj$Replicate, method = "absolute", 
                       supercell_membership = MC$membership)
    MC$integrated_clusters <- 
      supercell_assign(cluster = obj$integrated_clusters, method = "absolute", 
                       supercell_membership = MC$membership)
    MC$integrated_subclusters <- 
      supercell_assign(cluster = obj$integrated_subclusters, method = "absolute",
                       supercell_membership = MC$membership)
    MC$purity <- 
      supercell_purity(clusters = obj$integrated_clusters, 
                       method = "max_proportion", supercell_membership = MC$membership)
    
    # create the Seurat object version of this metacell object
    MC.seurat <- 
      supercell_2_Seurat( SC.GE = MC.ge, SC = MC, var.genes = MC$genes.use, N.comp = 20,
                          fields = c("integrated_clusters", 'Replicate',
                                     "integrated_subclusters", 'purity'))
    
    # transfer the other sample level metadata to the seurat object
    meta.rep2 = meta.rep[match(MC.seurat$Replicate, meta.rep$Replicate),]
    rownames(meta.rep2) = Cells(MC.seurat)
    MC.seurat = AddMetaData(MC.seurat, meta.rep2)
    
    # remove the metacells with low purity
    MC.seurat = subset(MC.seurat, purity > 0.5)
    
    saveRDS(MC.seurat, here(DATADIR, outs[i]))
    # }
}

###################################################################
## 2) merge together the supercells objects from different datasets

## combine the mouse and rat snRNA datasets from supercells
obj_list = lapply(here(DATADIR, outs[1:5]), readRDS)
obj_rodents = merge(obj_list[[1]], y = obj_list[-1])  
obj_rodents[['RNA']] = split(obj_rodents[['RNA']], obj_rodents$Project) 
obj_rodents = obj_rodents %>% FindVariableFeatures(verbose = F) %>% 
  ScaleData(verbose = F) %>% RunPCA(verbose = F)
obj_rodents = obj_rodents %>% 
  IntegrateLayers( method = HarmonyIntegration, orig.reduction = "pca", 
                   new.reduction = "integrated.rpca", verbose = F) %>% 
  RunUMAP(reduction = "integrated.rpca", verbose = F, dims = 1:30) %>% 
  JoinLayers()

obj_rodents %>%
  saveRDS(here('data/tidy_data/rdas/integrated_rodent_striatum.supercells.rds'))


## combine the human and macaque snRNA datasets from supercells
obj_list = lapply(here(DATADIR, outs[6:11]), readRDS)
obj_primate = merge(obj_list[[1]], y = obj_list[-1])
obj_primate[['RNA']] = split(obj_primate[['RNA']], obj_primate$Project) 
obj_primate = obj_primate %>% FindVariableFeatures(verbose = F) %>% 
  ScaleData(verbose = F) %>% RunPCA(verbose = F)
obj_primate = obj_primate %>% 
  IntegrateLayers( method = HarmonyIntegration, orig.reduction = "pca", 
                   new.reduction = "integrated.rpca", verbose = F) %>% 
  RunUMAP(reduction = "integrated.rpca", verbose = F, dims = 1:30) %>% 
  JoinLayers()

obj_primate %>% 
  saveRDS(here('data/tidy_data/rdas/integrated_primate_striatum.supercells.rds'))




