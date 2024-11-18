library(tidyverse)
library(Seurat)
library(data.table)
library(patchwork)
library(cluster)
library(future)
library(here)
library(AUCell)

# Enable parallelization
plan("multicore", workers = 12)
options(future.globals.maxSize = 100 * 1024^3, future.seed=TRUE)

DATADIR = 'data/tidy_data/'
PLOTDIR = 'figures/mouse_tshz1_pdyn/plots'
TABLDIR = 'figures/mouse_tshz1_pdyn/tables'
dir.create(PLOTDIR, recursive = T)
dir.create(TABLDIR, recursive = T)

##################################
## 0) prepare the metadata
proj1 = c('Gayden Kozel et al.','Siletti et al.', 'Tran et al.', 'Phan et al.', 
          'Li et al.', 'Corces et al.', 'Chiou et al.', 'He, Kleyman et al.', 
          'Phillips et al.', 'Savell et al.', 'Zeng et al.', 'Chen et al.',
          'Saunders et al.',  'Stanley et al.',  'Zu et al.')

clusters = c('D1', 'eSPN', 'D2', 'IC')
cols1 = setNames(RColorBrewer::brewer.pal(7,'Paired')[c(2,4,6,7)], clusters)
subclusters = c('D1-Matrix', 'D1-Striosome', 'D1-NUDAP',  'D1/D2H', 
                'D2-Matrix', 'D2-Striosome', 'IC')
cols2 = setNames(RColorBrewer::brewer.pal(7,'Paired'), subclusters)

#############################################################
## 1) load the single cell embedding across multiple species
multispecies_fn = 
  here('data/tidy_data/rdas/integrated_multispecies_striatum_ARC_P15.filtered.rds')
obj_mouse = readRDS(multispecies_fn) %>% subset(Project == 'Zeng et al.') %>% 
  subset(Region %in% 'NAcc') 

##################################################################
## 2) plot the integrated subclusters labels for the mouse NAcc
pdf(here(PLOTDIR, 'mouse_nac_integrated_subclusters.umap.pdf'), 
    height = 1.5, width = 1.5)
DimPlot(obj_mouse, group.by = 'integrated_subclusters', cols = cols2,
            raster = F, reduction = "umap.full", label = T, label.size = 2) & 
  FontSize(main = NA, x.title = 5, y.title = 5) & NoAxes() &
  theme_classic(base_size = 5) & 
  theme(legend.position = 'none', plot.title = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(),
        legend.title=element_text(size=4), legend.key.height = unit(.5, "lines"), 
        legend.text=element_text(size=4), legend.key.width  = unit(.2, "lines"), 
  )
dev.off()

##########################################################################
## 3) plot the dopamine receptors in conjunction with TSHZ1, PDYN, and OPRM1
pdf(here(PLOTDIR, 'mouse_tshz1_pdyn_oprm1_drd.umap.pdf'), height = 3, width = 2)
FeaturePlot(obj_mouse, features = c('DRD1', 'DRD2', 'TSHZ1', 'PDYN', 'OPRM1', 'PENK'),
            raster = F, reduction = "umap.full", cols = c('#cccccc50', '#DE2D26'), 
            ncol = 2,  min.cutoff = 'q5', max.cutoff = 'q90') & 
  FontSize(main = NA, x.title = 5, y.title = 5) & NoAxes() &
  theme_classic(base_size = 5) & 
  guides(fill = guide_colorbar(barwidth = 0.05, barheight = .7)) &
  theme(legend.position = 'right', 
        axis.ticks = element_blank(), axis.text = element_blank(),
        legend.title=element_text(size=4), legend.key.height = unit(.5, "lines"), 
        legend.text=element_text(size=4), legend.key.width  = unit(.2, "lines"), 
  )
dev.off()

