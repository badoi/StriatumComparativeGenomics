library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(data.table)
library(patchwork)
library(viridis)
library(cluster)
library(future)
library(here)

# Enable parallelization
plan("multicore", workers = 24)
options(future.globals.maxSize = 100 * 1024^3, future.seed=TRUE)

DATADIR = 'data/tidy_data/'
PLOTDIR = 'figures/multispecies_integration/plots'
dir.create(PLOTDIR, recursive = T)

#############################################################
## 1) load the single cell embedding across multiple species
save_fn = here(DATADIR,'rdas/integrated_mouse_rat_macaque_human_striatum.RNA.h5Seurat')
obj = LoadH5Seurat(save_fn)
Idents(obj) = 'Project'
DefaultAssay(obj) = 'RNA'

obj = obj

## from supplemental table s1 
tran_et_al_sex = c("M", 'M', "M", 'F', "M", 'M', "M", 'F')
names(tran_et_al_sex) = paste0('donor', 1:8)

hekleyman_et_al_sex = c('Monkey_F' = 'F', 'Monkey_P' = 'M')

obj@meta.data$Sex = 
  ifelse(obj$Project == 'Phan_et_al', obj$Sex,  
         ifelse(obj$Project == 'Stanley_et_al', 'M', 
                c(tran_et_al_sex, hekleyman_et_al_sex)[obj$Replicate]))

# ################################################################################
# # 3) Violin plot - Visualize single cell expression distributions in each cluster
# obj_tran = subset(obj, integrated_clusters != 'Other')
# obj_tran = subset(obj_tran, Project == 'Tran_et_al')
# obj_tran$integrated_clusters2 = 
#   ifelse(obj_tran$integrated_clusters == 'D1/2-Hybrid',  
#          'D1/2H', obj_tran$integrated_clusters)
# 
# features = c('DRD3')
# pdf(here::here(PLOTDIR, 'sex_specific_human_DRD3.vln.pdf'), width = 1.5, height = 2)
# VlnPlot(obj_tran, features = features, assay = 'RNA', group.by = 'integrated_clusters2', 
#         ncol = 1, pt.size = F, split.by = 'Sex', split.plot = TRUE) & 
#   facet_grid(~factor(obj_tran$Project, levels)) & 
#   theme(axis.title.x = element_blank()) &
#   FontSize(main = 7, x.title = 5, y.title = 5, x.text = 5, y.text = 5) & 
#   theme(strip.text.x = element_text(size = 5),
#         legend.position = 'bottom', 
#         axis.title.y = element_blank(), axis.text.x = element_text(angle = 0), 
#         plot.title = element_blank()
#         )
# dev.off()


