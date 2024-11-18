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
cols_groups = setNames(RColorBrewer::brewer.pal(7,'Paired'), subclusters)

mytheme = FontSize(main = NA, x.title = 5, y.title = 5) & NoAxes() &
  theme_classic(base_size = 5) & 
  theme(legend.position = 'none', 
        legend.title=element_text(size=4), legend.key.height = unit(.5, "lines"), 
        legend.text=element_text(size=4), legend.key.width  = unit(.2, "lines"))

#############################################################
## 1) create and load the pseudo bulk object for the analysis
supercells_fn = here('data/tidy_data/Mouse_Zeng/WMB-10Xv3-STR.msn.seurat.supercells.rds')
obj_mouse = readRDS(supercells_fn) %>% subset(Region %in% 'NAcc') 

## add the Pdyn/Tshz1 expression thresholds
thresh1 = 1.5
keep_groups = c('Tshz1+, Pdyn-', 'Tshz1-, Pdyn+', 
                'Tshz1-, Pdyn-', 'Tshz1+, Pdyn+')
cols_groups = c('Tshz1+, Pdyn-' = '#b856d7', 'Tshz1-, Pdyn+' = '#8babd3', 
                'Tshz1-, Pdyn-' = 'grey', 'Tshz1+, Pdyn+' = 'black')

df = LayerData(obj_mouse, 'data')[c('Tshz1', 'Pdyn'),] %>% as.matrix() %>% 
  t() %>% as.data.frame() %>% bind_cols(x = obj_mouse[[]]) %>% 
  mutate(group = 
           case_when(Tshz1 > thresh1 & Pdyn < thresh1 ~ 'Tshz1+, Pdyn-',
                     Pdyn > thresh1 & Tshz1 < thresh1~ 'Tshz1-, Pdyn+',
                     Pdyn < thresh1 & Tshz1 < thresh1 ~ 'Tshz1-, Pdyn-', 
                     Pdyn > thresh1 & Tshz1 > thresh1 ~ 'Tshz1+, Pdyn+'),
         group = factor(group, levels = keep_groups))

obj_mouse$group = df$group
table(df$group)/nrow(df) * 100 %>% signif(2)
table(obj_mouse$integrated_subclusters, obj_mouse$group)

# violin plot of the expression of serotonin receptors
genes_mm = grep(Features(obj_mouse), pattern = 'Htr[0-9]', value = T) %>% sort()

pdf(here(PLOTDIR, 'mouse_nacc_pdyn_tshz1_htr.vln.pdf'), 
    height = 6, width = 6)
VlnPlot(obj_mouse, genes_mm, group.by = 'group', cols = cols_groups, 
        raster = T, pt.size = 0, same.y.lims = T) & mytheme & coord_flip() & 
  scale_x_discrete(limits = rev)
dev.off()


#############################################################
## 2) load the single cell embedding across multiple species
obj_human = here('data/tidy_data/rdas/integrated_primate_striatum.supercells.rds') %>% 
  readRDS() %>% subset(Region %in% 'NAcc') %>% subset(Species == 'Human') 

thresh1 = 1; thresh2 = .5  
keep_groups = c('TSHZ1+, PDYN-', 'TSHZ1-, PDYN+', 'TSHZ1-, PDYN-', 'TSHZ1+, PDYN+')
cols_groups = c('TSHZ1+, PDYN-' = '#b856d7', 'TSHZ1-, PDYN+' = '#8babd3', 
                'TSHZ1-, PDYN-' = 'grey', 'TSHZ1+, PDYN+' = 'black')

df2 = LayerData(obj_human, 'data')[c('TSHZ1', 'PDYN'),] %>% as.matrix() %>% 
  t() %>% as.data.frame() %>% bind_cols(x = obj_human[[]]) %>% 
  mutate(group = 
           case_when(TSHZ1 > thresh1 & PDYN < thresh2 ~ 'TSHZ1+, PDYN-',
                     PDYN > thresh2 & TSHZ1 < thresh1 ~ 'TSHZ1-, PDYN+',
                     PDYN < thresh2 & TSHZ1 < thresh1 ~ 'TSHZ1-, PDYN-', 
                     PDYN > thresh2 & TSHZ1 > thresh1 ~ 'TSHZ1+, PDYN+'),
         group = factor(group, levels = keep_groups))

obj_human$group = df2$group
table(df2$group)/nrow(df2) * 100 %>% signif(2)
table(obj_human$integrated_subclusters, obj_human$group)

# violin plot of the expression of serotonin receptors
genes_hs = grep(Features(obj_human), pattern = 'HTR[0-9]', value = T) %>% sort()

pdf(here(PLOTDIR, 'human_nacc_pdyn_tshz1_htr.vln.pdf'), 
    height = 8, width = 6)
VlnPlot(obj_human, genes_hs, group.by = 'group', cols = cols_groups, 
        raster = T, pt.size = 0, same.y.lims = T) & mytheme & coord_flip() & 
  scale_x_discrete(limits = rev)
dev.off()

