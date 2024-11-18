library(tidyverse)
library(Seurat)
library(data.table)
library(limma)
library(edgeR)
library(future)
library(here)

# Enable parallelization
plan("multicore", workers = 12)
options(future.globals.maxSize = 100 * 1024^3, future.seed=TRUE)

DATADIR = 'data/tidy_data/'
PLOTDIR = 'figures/mouse_tshz1_pdyn/plots'
dir.create(PLOTDIR, recursive = T)

TABLDIR = 'figures/mouse_tshz1_pdyn/tables'
dir.create(TABLDIR, recursive = T)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

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

mytheme = FontSize(main = NA, x.title = 5, y.title = 5) & NoAxes() &
  theme_classic(base_size = 5) & 
  theme(legend.position = 'none', 
        legend.title=element_text(size=4), legend.key.height = unit(.5, "lines"), 
        legend.text=element_text(size=4), legend.key.width  = unit(.2, "lines"))

#############################################################
## 1) create and load the pseudo bulk object for the analysis
supercells_fn = here('data/tidy_data/Mouse_Zeng/WMB-10Xv3-STR.msn.seurat.supercells.rds')
obj_mouse = readRDS(supercells_fn) %>% subset(Region %in% 'NAcc') 

# violin plot of the expression of serotonin receptors
genes_mm = grep(Features(obj_mouse), pattern = 'Htr[0-9]', value = T) %>% sort()

pdf(here(PLOTDIR, 'mouse_nacc_htr.vln.pdf'), 
    height = 6, width = 6)
VlnPlot(obj_mouse, genes_mm, group.by = 'integrated_subclusters', cols = cols2, 
        raster = T, pt.size = 0, same.y.lims = T) & mytheme & coord_flip() & 
  scale_x_discrete(limits = rev)
dev.off()


#############################################################
## 2) load the single cell embedding across multiple species
obj_human = here('data/tidy_data/rdas/integrated_primate_striatum.supercells.rds') %>% 
  readRDS() %>% subset(Region %in% 'NAcc') %>% subset(Species == 'Human') 

# violin plot of the expression of serotonin receptors
genes_hs = grep(Features(obj_human), pattern = 'HTR[0-9]', value = T) %>% sort()

pdf(here(PLOTDIR, 'human_nacc_htr.vln.pdf'), 
    height = 8, width = 6)
VlnPlot(obj_human, genes_hs, group.by = 'integrated_subclusters', cols = cols2, 
        raster = T, pt.size = 0, same.y.lims = T) & mytheme & coord_flip() & 
  scale_x_discrete(limits = rev)
dev.off()

